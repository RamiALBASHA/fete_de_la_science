from datetime import datetime
from json import dump
from json import load as load_json
from pathlib import Path
from pickle import load as load_pickle

from hydroshoot import architecture, display, io, initialisation, model, constants, exchange
from k3d.plot import Plot
from matplotlib import pyplot, patches
from matplotlib.dates import DateFormatter
from oawidgets.plantgl import PlantGL
from openalea.mtg import traversal, mtg
from openalea.plantgl.all import Scene, Viewer
from openalea.plantgl.all import surface as surf
from pandas import read_csv, to_datetime

pyplot.ioff()

CONV_CO2 = constants.co2_molar_mass * 1.e-6 * 3600.  # umol/s to g/h

MAP_NAMES = dict(
    vsp='Espalier bas',
    lyre='Lyre',
    gdc='Rideau simple',
    cordon='Cordon libre')

PATH_STATIC = Path(__file__).parent / 'data'


def build_grapevine_mtg(path_digit_file: Path, is_cordon_preferential_orientation: bool = False) -> mtg.MTG:
    g = architecture.vine_mtg(file_path=path_digit_file)
    cordon_vector = architecture.cordon_vector(g=g)[1] if is_cordon_preferential_orientation else None
    # Local Coordinates Correction
    for vid in traversal.iter_mtg2(g, g.root):
        # print(v)
        architecture.vine_phyto_modular(g, vid)
        # if is_high_leaf_area:
        #     architecture.vine_axeII(g, v, pruning_type='avg_field_model', N_max=6, insert_angle=90, N_max_order=6,
        #                             insert_angle_CI=0, a_L=80)
        architecture.vine_petiole(g, vid, pet_ins=90., phyllo_angle=180, phyllo_angle_cv=0)
        architecture.vine_leaf(g, vid, leaf_inc=-45., lim_max=12.5, rand_rot_angle=30,
                               lim_min=5., order_lim_max=6., max_order=55, cordon_vector=cordon_vector)
        architecture.vine_mtg_properties(g, vid)
        architecture.vine_mtg_geometry(g, vid)
        architecture.vine_transform(g, vid)

    return g


def generate_grapevine_mtgs(path_dir_project: Path):
    path_output = path_dir_project / 'preprocessed_inputs'
    for is_high_leaf_area in (False, True):
        g = build_grapevine_mtg(
            path_digit_file=path_dir_project / f'cordon_bilateral_{"high" if is_high_leaf_area else "low"}.csv')
        virtual_scene = display.visu(g=g, def_elmnt_color_dict=True, scene=Scene(), view_result=True)
        total_leaf_area = calc_total_leaf_area(g=g)
        preprocess_inputs(
            grapevine_mtg=g, path_project_dir=path_dir_project, psi_soil=-1, gdd_since_budbreak=1000.,
            scene=virtual_scene, name_index=f'{total_leaf_area:.2f}')
        architecture.save_mtg(
            g=g, scene=virtual_scene, file_path=path_output, filename=f'mtg{total_leaf_area:.2f}.pckl')
    pass


def generate_grapevine_mtgs2(path_dir_project: Path):
    path_output = path_dir_project / 'preprocessed_inputs'
    for training_sys in ('cordon', 'lyre', 'gdc', 'vsp'):
        g = build_grapevine_mtg(path_digit_file=path_dir_project / f'digit_{training_sys}.csv')
        virtual_scene = display.visu(g=g, def_elmnt_color_dict=True, scene=Scene(), view_result=True)
        preprocess_inputs(
            grapevine_mtg=g, path_project_dir=path_dir_project, psi_soil=-1, gdd_since_budbreak=1000.,
            scene=virtual_scene, name_index=training_sys)
        architecture.save_mtg(
            g=g, scene=virtual_scene, file_path=path_output, filename=f'mtg{training_sys}.pckl')
    pass


def calc_total_leaf_area(g: mtg.MTG) -> float:
    total_leaf_area = 0
    for vid in g.VtxList(Scale=3):
        n = g.node(vid)
        if n.label.startswith('L'):
            total_leaf_area += surf(n.geometry)
    return total_leaf_area * 1.e-4


def preprocess_inputs(grapevine_mtg: mtg.MTG, path_project_dir: Path, psi_soil: float, gdd_since_budbreak: float,
                      scene: Scene, name_index: str = None):
    inputs = io.HydroShootInputs(g=grapevine_mtg, path_project=path_project_dir, scene=scene, psi_soil=psi_soil,
                                 gdd_since_budbreak=gdd_since_budbreak)
    io.verify_inputs(g=grapevine_mtg, inputs=inputs)
    grapevine_mtg = initialisation.init_model(g=grapevine_mtg, inputs=inputs)

    static_data = {'form_factors': {s: grapevine_mtg.property(s) for s in ('ff_sky', 'ff_leaves', 'ff_soil')}}
    static_data.update({'Na': grapevine_mtg.property('Na')})

    dynamic_data = {}
    inputs_hourly = io.HydroShootHourlyInputs(psi_soil=inputs.psi_soil_forced, sun2scene=inputs.sun2scene)
    for date_sim in inputs.params.simulation.date_range:
        print(date_sim)
        inputs_hourly.update(
            g=grapevine_mtg, date_sim=date_sim, hourly_weather=inputs.weather[inputs.weather.index == date_sim],
            psi_pd=inputs.psi_pd, params=inputs.params)

        grapevine_mtg, diffuse_to_total_irradiance_ratio = initialisation.init_hourly(
            g=grapevine_mtg, inputs_hourly=inputs_hourly, leaf_ppfd=inputs.leaf_ppfd,
            params=inputs.params)

        dynamic_data.update({grapevine_mtg.date: {
            'diffuse_to_total_irradiance_ratio': diffuse_to_total_irradiance_ratio,
            'Ei': grapevine_mtg.property('Ei'),
            'Eabs': grapevine_mtg.property('Eabs')}})

    path_preprocessed_inputs = path_project_dir / 'preprocessed_inputs'
    path_preprocessed_inputs.mkdir(exist_ok=True)
    with open(path_preprocessed_inputs / f'static{name_index}.json', mode='w') as f_prop:
        dump(static_data, f_prop)
    pass
    with open(path_preprocessed_inputs / f'dynamic{name_index}.json', mode='w') as f_prop:
        dump(dynamic_data, f_prop)
    pass


def run_hydroshoot(path_dir_preprocessed_inputs: Path, path_project: Path, path_output: Path, training_sys: str,
                   is_show: bool = True, snap_shot_path: Path = False):
    path_mtg_file = path_dir_preprocessed_inputs / f'mtg{training_sys}.pckl'
    with open(path_mtg_file, mode='rb') as f:
        g, _ = load_pickle(f)

    # eval('%gui qt5')
    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_mtg_geometry(g, v)
        architecture.vine_transform(g, v)

    virtual_scene = display.visu(g=g, def_elmnt_color_dict=True, scene=Scene(), view_result=False)
    if is_show:
        Viewer.display(virtual_scene)
    if snap_shot_path is not None:
        Viewer.saveSnapshot(str(snap_shot_path))

    with open(path_dir_preprocessed_inputs / f'static{training_sys}.json') as f:
        static_inputs = load_json(f)
    with open(path_dir_preprocessed_inputs / f'dynamic{training_sys}.json') as f:
        dynamic_inputs = load_json(f)

    summary_output = model.run(g=g, wd=path_project, scene=virtual_scene,
                               psi_soil=-0.5, gdd_since_budbreak=1000.,
                               form_factors=static_inputs['form_factors'],
                               leaf_nitrogen=static_inputs['Na'],
                               leaf_ppfd=dynamic_inputs,
                               path_output=path_output)

    return g, summary_output


def display_whole_plant(path_output_low: Path, path_output_high: Path, path_weather: Path) -> pyplot.Figure:
    df_output_low = read_csv(path_output_low, sep=';', decimal='.', index_col='Unnamed: 0')
    df_output_high = read_csv(path_output_high, sep=';', decimal='.', index_col='Unnamed: 0')
    for df in (df_output_low, df_output_high):
        df.index = to_datetime(df.index)

    df_weather = read_csv(path_weather, sep=';', decimal='.', index_col='time')
    df_weather.index = to_datetime(df_weather.index)
    df_weather = df_weather.loc[df_output_low.index, :]

    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharex='all', figsize=(8, 4))

    axs[0, 0].plot(df_output_low.index, df_weather['Rg'], label='Incident', c='k')
    axs[0, 0].plot(df_output_low.index, df_output_low['Rg'], label='Absorbé (canopée clairsemée)', c='orange')
    axs[0, 0].plot(df_output_high.index, df_output_high['Rg'], label='Absorbé (canopée dense)', c='red')
    axs[0, 0].set(ylabel='Rayonnement solaire\n' + '($\mathregular{W\/m^{-2}_{sol}}$)')
    axs[0, 0].legend(loc='center left', fontsize=8)

    axs[1, 0].plot(df_output_low.index, df_weather['Tac'], label='air', c='k')
    axs[1, 0].plot(df_output_low.index, df_output_low['Tleaf'], label='Canopée (clairsemée)', c='orange')
    axs[1, 0].plot(df_output_high.index, df_output_high['Tleaf'], label='Canopée (dense)', c='red')
    axs[1, 0].set(ylabel='Température\n(°C)', xlabel='heure')
    axs[1, 0].legend(loc='lower right', fontsize=8)

    axs[0, 1].plot(df_output_low.index, df_output_low['An'] * CONV_CO2, c='orange', label='(clairsemée)')
    axs[0, 1].plot(df_output_high.index, df_output_high['An'] * CONV_CO2, c='red', label='dense')
    axs[0, 1].set(ylabel='\n'.join(('Photosynthèse', r'($\mathregular{g\/plant^{-1}}$)')))
    axs[0, 1].legend(loc='upper left', fontsize=8)

    axs[1, 1].plot(df_output_low.index, df_output_low['E'], c='orange', label='(clairsemée)')
    axs[1, 1].plot(df_output_high.index, df_output_high['E'], c='red', label='dense')
    axs[1, 1].set(ylabel='\n'.join(('Transpiration', r'($\mathregular{g\/plant^{-1}}$)')), xlabel='heure')
    axs[1, 1].legend(loc='upper left', fontsize=8)

    axs[-1, -1].xaxis.set_major_formatter(DateFormatter('%H'))
    fig.tight_layout()
    return fig


def plot_whole_plant_gas_exchange(path_ref: Path, path_user: Path, training_system_1: str, training_system_2: str,
                                  training_system_user: str, path_weather: Path) -> pyplot.Figure:
    path_output_1 = path_ref / training_system_1 / f'time_series.csv'
    path_output_2 = path_ref / training_system_2 / f'time_series.csv'
    df_output_ref_1 = read_csv(path_output_1, sep=';', decimal='.', index_col='Unnamed: 0')
    df_output_ref_2 = read_csv(path_output_2, sep=';', decimal='.', index_col='Unnamed: 0')
    for df in (df_output_ref_1, df_output_ref_2):
        df.index = to_datetime(df.index)

    df_weather = read_csv(path_weather, sep=';', decimal='.', index_col='time')
    df_weather.index = to_datetime(df_weather.index)
    df_weather = df_weather.loc[df_output_ref_1.index, :]

    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharex='all')

    axs[0, 0].plot(df_weather.index, df_weather['Rg'], label='Incident', c='k')
    axs[0, 0].plot(df_output_ref_1.index, df_output_ref_1['Rg'], c='orange', label=MAP_NAMES[training_system_1])
    axs[0, 0].plot(df_output_ref_2.index, df_output_ref_2['Rg'], c='red', label=MAP_NAMES[training_system_2])
    axs[0, 0].set(ylabel='Rayonnement solaire\n' + '($\mathregular{W\/m^{-2}_{sol}}$)')

    axs[0, 1].plot(df_output_ref_1.index, df_output_ref_1['An'] * CONV_CO2, c='orange',
                   label=MAP_NAMES[training_system_1])
    axs[0, 1].plot(df_output_ref_2.index, df_output_ref_2['An'] * CONV_CO2, c='red', label=MAP_NAMES[training_system_2])
    axs[0, 1].set(ylabel='\n'.join(('Photosynthèse nette', r'($\mathregular{g\/plant^{-1}}$)')))

    axs[1, 1].plot(df_output_ref_1.index, df_output_ref_1['E'], c='orange', label=MAP_NAMES[training_system_1])
    axs[1, 1].plot(df_output_ref_2.index, df_output_ref_2['E'], c='red', label=MAP_NAMES[training_system_2])
    axs[1, 1].set(ylabel='\n'.join(('Transpiration', r'($\mathregular{g\/plant^{-1}}$)')))

    axs[1, 0].plot(df_output_ref_1.index, df_weather['Tac'], label='air', c='k')
    axs[1, 0].plot(df_output_ref_1.index, df_output_ref_1['Tleaf'], label=MAP_NAMES[training_system_1], c='orange')
    axs[1, 0].plot(df_output_ref_2.index, df_output_ref_2['Tleaf'], label=MAP_NAMES[training_system_2], c='red')
    axs[1, 0].set(ylabel='Température\n(°C)', xlabel='heure')

    if training_system_user != "none":
        path_output_user = path_user / f'time_series.csv'
        df_output_user = read_csv(path_output_user, sep=';', decimal='.', index_col='Unnamed: 0')
        df_output_user.index = to_datetime(df_output_user.index)
        kwargs = dict(label=MAP_NAMES[training_system_user], c='blue')
        axs[0, 0].plot(df_output_user.index, df_output_user['Rg'], **kwargs)
        axs[0, 1].plot(df_output_user.index, df_output_user['An'] * CONV_CO2, **kwargs)
        axs[1, 1].plot(df_output_user.index, df_output_user['E'], **kwargs)
        axs[1, 0].plot(df_output_user.index, df_output_user['Tleaf'], **kwargs)
        pass

    axs[0, 0].legend()
    axs[1, 0].legend(loc='lower right', fontsize=8)

    for ax in axs[-1, :]:
        ax.xaxis.set_major_formatter(DateFormatter('%H'))
    fig.tight_layout()
    return fig


def plot_water_use_efficiency(path_ref: Path, path_user: Path, training_system_1: str, training_system_2: str,
                              training_system_user: str) -> pyplot.Figure:
    path_output_1 = path_ref / training_system_1 / f'time_series.csv'
    path_output_2 = path_ref / training_system_2 / f'time_series.csv'
    df_ref_1 = read_csv(path_output_1, sep=';', decimal='.', index_col='Unnamed: 0')
    df_ref_1.index = to_datetime(df_ref_1.index)

    df_ref_2 = read_csv(path_output_2, sep=';', decimal='.', index_col='Unnamed: 0')
    df_ref_2.index = to_datetime(df_ref_2.index)

    xticklebels = [MAP_NAMES[s] for s in (training_system_1, training_system_2)]
    fig, ax = pyplot.subplots()
    for x, df, c in ((0, df_ref_1, 'orange'), (1, df_ref_2, 'red')):
        ax.bar(x, df.loc[df_ref_2.index, 'An'].sum() * CONV_CO2 / df['E'].sum(), facecolor=c)

    if training_system_user != 'none':
        df_user = read_csv(path_user / 'time_series.csv', sep=';', decimal='.', index_col='Unnamed: 0')
        df_user.index = to_datetime(df_user.index)
        ax.bar(2, df_user.loc[df_user.index, 'An'].sum() * CONV_CO2 / df_user['E'].sum(),
               facecolor='blue')
        xticklebels.append(MAP_NAMES[training_system_user])

    ax.set_xticks(range(len(xticklebels)))
    ax.set_xticklabels(xticklebels)

    ax.set(ylabel='\n'.join(["Efficience de l'utilisation de l'eau", "$\mathregular{g_{CO_2}\/g^{-1}_{H_2O}}$"]))
    fig.tight_layout()

    return fig


def display_mtg_properties(path_output_low: Path, path_output_high: Path):
    fig, axs = pyplot.subplots(nrows=4, ncols=1, sharex='all')
    irradiance = {}
    photosynthesis = {}
    transpiration = {}
    temperature = {}
    nb_leaves = {}
    for s, path_dir in (('low', path_output_low), ('high', path_output_high)):
        pckl_files = [f for f in path_dir.iterdir() if f.name.endswith('.pckl')]
        for pckl_file in pckl_files:
            with open(pckl_file, mode='rb') as f:
                g, _ = load_pickle(f)
            date_sim = datetime.strptime(g.date, "%Y%m%d%H%M%S")
            if date_sim not in irradiance:
                irradiance[date_sim] = {}
            irradiance[date_sim][s] = list(g.property('Rg').values())
            if date_sim not in photosynthesis:
                photosynthesis[date_sim] = {}
            photosynthesis[date_sim][s] = [v * CONV_CO2 for v in g.property('FluxC').values()]
            if date_sim not in transpiration:
                transpiration[date_sim] = {}
            transpiration[date_sim][s] = [v * 3600. for v in g.property('Flux').values()]
            if date_sim not in temperature:
                temperature[date_sim] = {}
            temperature[date_sim][s] = list(g.property('Tlc').values())
            nb_leaves.update({s: len(g.property('An'))})

    bplots = {s: {'low': [], 'high': []} for s in ('irradiance', 'photosynthesis', 'transpiration', 'temperature')}
    for ax, (var_name, var_value) in zip(axs, (('irradiance', irradiance), ('photosynthesis', photosynthesis),
                                               ('transpiration', transpiration), ('temperature', temperature))):
        for sim_datetime in var_value.keys():
            for leaf_area in ('low', 'high'):
                x_position = sim_datetime.hour + (-0.1 if leaf_area == 'low' else +0.1)
                bplot = ax.boxplot(x=var_value[sim_datetime][leaf_area], patch_artist=True, vert=True,
                                   positions=[x_position], widths=0.1, sym='')
                bplot['boxes'][0].set_facecolor('orange' if leaf_area == 'low' else 'red')
                bplot['boxes'][0].set_edgecolor('orange' if leaf_area == 'low' else 'red')
                bplots[var_name][leaf_area].append(bplot)
    axs[0].set(ylabel='\n'.join(('Rayonnement absorbé', '($\mathregular{W\/m^{-2}_{feuille}}$)')))
    axs[1].set(ylabel='\n'.join(('Photosynthèse nette', r'($\mathregular{g\/{feuille}^{-1}}$)')))
    axs[2].set(ylabel='\n'.join(('Transpiration', r'($\mathregular{g\/{feuille}^{-1}}$)')))
    axs[3].set(ylabel='\n'.join(('Température', '(°C)')), xlabel='heure')

    axs[0].text(0.025, 0.7, f'Nb feuilles:\n clairsemé: {nb_leaves["low"]}\n dense: {nb_leaves["high"]}',
                transform=axs[0].transAxes)

    handles_labels = [patches.Patch(facecolor='orange' if s == 'clairsemé' else 'red', label=s)
                      for s in ('clairsemé', 'dense')]
    axs[0].legend(handles=handles_labels)
    pass


def display_mtg_property(path_ref: Path, training_system_1: str, training_system_2: str, path_user: Path,
                         training_system_user: str, hour: int, mtg_property: str) -> pyplot.Figure:
    path_ref_1 = path_ref / training_system_1 / f'mtg20120801{hour:02d}0000.pckl'
    path_ref_2 = path_ref / training_system_2 / f'mtg20120801{hour:02d}0000.pckl'
    with open(path_ref_1, mode='rb') as f:
        g_ref_1, _ = load_pickle(f)
    with open(path_ref_2, mode='rb') as f:
        g_ref_2, _ = load_pickle(f)

    n_cols = 2
    if training_system_user != 'none':
        n_cols = 3
        path_user = path_user / f'mtg20120801{hour:02d}0000.pckl'
        with open(path_user, mode='rb') as f:
            g_user, _ = load_pickle(f)

    if mtg_property == 'psi_head':
        prop2 = None
        cmap = None
    else:
        prop2 = 'Eabs'
        cmap = 'inferno'

    fig, axs = pyplot.subplots(ncols=n_cols, sharey='all', sharex='all')
    axs[0] = display.property_map(g_ref_1, prop=mtg_property, ax=axs[0], prop2=prop2, colormap=cmap)
    axs[0].set_title(MAP_NAMES[training_system_1])
    axs[1] = display.property_map(g_ref_2, prop=mtg_property, ax=axs[1], prop2=prop2, colormap=cmap)
    axs[1].set_title(MAP_NAMES[training_system_2])
    if training_system_user != 'none':
        axs[2] = display.property_map(g_user, prop=mtg_property, ax=axs[2], prop2=prop2, colormap=cmap)
        axs[2].set_title(MAP_NAMES[training_system_user])
    return fig


def display_water_use_efficiency(path_output_low: Path, path_output_high: Path) -> pyplot.Figure:
    df_output_low = read_csv(path_output_low, sep=';', decimal='.', index_col='Unnamed: 0')
    df_output_high = read_csv(path_output_high, sep=';', decimal='.', index_col='Unnamed: 0')

    fig, ax = pyplot.subplots()
    for s, df, c in (('clairsemé', df_output_low, 'orange'), ('dense', df_output_high, 'red')):
        df.index = to_datetime(df.index)
        ax.bar(s, df['An'].sum() * CONV_CO2 / df['E'].sum(), facecolor=c)
    ax.set(ylabel='\n'.join(["Efficience de l'utilisation de l'eau", "$\mathregular{g_{CO_2}\/g^{-1}_{H_2O}}$"]))
    fig.tight_layout()

    return fig


def update_params(
        path_default_params: Path, path_user_params: Path, stomatal_conductance_min: float, stomatal_psi50: float,
        steepness: float):
    with open(path_default_params, mode='r') as f:
        params = load_json(f)
    params['simulation']['sdate'] = "2012-08-01 09:00:00"
    params['simulation']['edate'] = "2012-08-01 15:00:00"
    params['exchange']['par_gs']['g0'] = stomatal_conductance_min * 1.e-3
    params['exchange']['par_gs']['psi0'] = stomatal_psi50
    params['exchange']['par_gs']['n'] = steepness
    with open(path_user_params, mode='w') as f:
        dump(params, f, indent=2)
    pass


def plot_stomatal_reduction_coef(psi50: float, gs0: float, n: float) -> pyplot.Figure:
    gs0 /= 0.4
    fig, ax = pyplot.subplots()
    psi_leaf = [v / 1000 for v in range(-4000, 0, 10)]
    gs_ref = [exchange.fvpd_3(model='misson', psi=v, vpd=None, psi_crit=-1, steepness_tuzet=4, m0=1) + gs0
              for v in psi_leaf]
    gs = [exchange.fvpd_3(model='misson', psi=v, vpd=None, psi_crit=psi50, steepness_tuzet=n, m0=1) + gs0
          for v in psi_leaf]
    ax.plot(psi_leaf, gs_ref, label='ref', color='grey')
    ax.plot(psi_leaf, gs, label='usr', color='blue')
    ax.vlines(x=-1, ymin=0, ymax=gs0 + 0.5, color='grey', linestyles='--')
    ax.vlines(x=psi50, ymin=0, ymax=gs0 + 0.5, color='blue', linestyles='--')
    ax.hlines(y=0.02 / 0.4, xmin=min(psi_leaf), xmax=max(psi_leaf), color='grey', linestyles='--')
    ax.hlines(y=gs0, xmin=min(psi_leaf), xmax=max(psi_leaf), color='blue', linestyles='--')
    ax.set_yticks([0, gs0, 0.5 + gs0, max(gs)])
    ax.set_yticklabels(['0', 'cond. res.', 'cond. @ 50%', 'cond. max.'])
    ax.set_xlabel('Pression (MPa)')
    fig.tight_layout()
    return fig


def show_3d(canopy_name: str) -> Plot:
    g = build_grapevine_mtg(path_digit_file=PATH_STATIC / f'digit_{canopy_name}.csv')
    virtual_scene = display.visu(g=g, def_elmnt_color_dict=True, scene=Scene(), view_result=False)
    return PlantGL(virtual_scene)


if __name__ == '__main__':
    path_data = Path(__file__).parent / 'data'
    # generate_grapevine_mtgs(path_dir_project=path_data)
    generate_grapevine_mtgs2(path_dir_project=path_data)
    # %gui qt5

    for training_system in ('cordon', 'lyre', 'gdc', 'vsp'):
        vine_mtg, output = run_hydroshoot(
            path_project=path_data,
            path_dir_preprocessed_inputs=path_data / 'preprocessed_inputs',
            training_sys=training_system,
            is_show=False,
            path_output=path_data / f'output_ref/{training_system}/time_series.csv')
    #
    # fig1 = display_whole_plant(
    #     path_output_low=path_data / 'output_low/time_series.csv',
    #     path_output_high=path_data / 'output_high/time_series.csv',
    #     path_weather=path_data / 'weather.csv')
    #
    # display_mtg_properties(
    #     path_output_low=path_data / 'output_low',
    #     path_output_high=path_data / 'output_high')
    #
    # display_water_use_efficiency(
    #     path_output_low=path_data / 'output_low/time_series.csv',
    #     path_output_high=path_data / 'output_high/time_series.csv')
