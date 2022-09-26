from datetime import datetime
from json import dump
from json import load as load_json
from pathlib import Path
from pickle import load as load_pickle

from hydroshoot import architecture, display, io, initialisation, model, constants
from matplotlib import pyplot, patches
from matplotlib.dates import DateFormatter
from openalea.mtg import traversal, mtg
from openalea.plantgl.all import Scene
from openalea.plantgl.all import surface as surf
from pandas import read_csv, to_datetime

CONV_CO2 = constants.co2_molar_mass * 1.e-6 * 3600.  # umol/s to g/h


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


def run_hydroshoot(path_dir_preprocessed_inputs: Path, path_output: Path, is_low_leaf_area: bool, is_show: bool = True):
    leaf_area = 1.54 if is_low_leaf_area else 3.78
    path_mtg_file = path_dir_preprocessed_inputs / f'mtg{leaf_area}.pckl'
    with open(path_mtg_file, mode='rb') as f:
        g, _ = load_pickle(f)

    if is_show:
        # eval('%gui qt5')
        for v in traversal.iter_mtg2(g, g.root):
            architecture.vine_mtg_geometry(g, v)
            architecture.vine_transform(g, v)
        virtual_scene = display.visu(g=g, def_elmnt_color_dict=True, scene=Scene(), view_result=True)
    else:
        virtual_scene = None

    with open(path_dir_preprocessed_inputs / f'static{leaf_area}.json') as f:
        static_inputs = load_json(f)
    with open(path_dir_preprocessed_inputs / f'dynamic{leaf_area}.json') as f:
        dynamic_inputs = load_json(f)

    summary_output = model.run(g=g, wd=path_data, scene=virtual_scene,
                               psi_soil=-0.5, gdd_since_budbreak=1000.,
                               form_factors=static_inputs['form_factors'],
                               leaf_nitrogen=static_inputs['Na'],
                               leaf_ppfd=dynamic_inputs,
                               path_output=path_output)

    return g, summary_output


def display_whole_plant(path_output_low: Path, path_output_high: Path, path_weather: Path):
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
    axs[0, 0].plot(df_output_low.index, df_output_high['Rg'], label='Absorbé (canopée dense)', c='red')
    axs[0, 0].set(ylabel='Rayonnement solaire\n' + '($\mathregular{W\/m^{-2}_{sol}}$)')
    axs[0, 0].legend(loc='center left', fontsize=8)

    axs[1, 0].plot(df_output_low.index, df_weather['Tac'], label='air', c='k')
    axs[1, 0].plot(df_output_low.index, df_output_low['Tleaf'], label='Canopée (clairsemée)', c='orange')
    axs[1, 0].plot(df_output_low.index, df_output_high['Tleaf'], label='Canopée (dense)', c='red')
    axs[1, 0].set(ylabel='Température\n(°C)', xlabel='heure')
    axs[1, 0].legend(loc='lower right', fontsize=8)

    axs[0, 1].plot(df_output_low.index, df_output_low['An'] * CONV_CO2, c='orange', label='(clairsemée)')
    axs[0, 1].plot(df_output_low.index, df_output_high['An'] * CONV_CO2, c='red', label='dense')
    axs[0, 1].set(ylabel='\n'.join(('Photosynthèse nette', r'($\mathregular{g\/plant^{-1}}$)')))
    axs[0, 1].legend(loc='upper left', fontsize=8)

    axs[1, 1].plot(df_output_low.index, df_output_low['E'], c='orange', label='(clairsemée)')
    axs[1, 1].plot(df_output_low.index, df_output_high['E'], c='red', label='dense')
    axs[1, 1].set(ylabel='\n'.join(('Transpiration', r'($\mathregular{g\/plant^{-1}}$)')), xlabel='heure')
    axs[1, 1].legend(loc='upper left', fontsize=8)

    axs[-1, -1].xaxis.set_major_formatter(DateFormatter('%H'))
    fig.tight_layout()
    pass


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


def display_water_use_efficiency(path_output_low: Path, path_output_high: Path):
    df_output_low = read_csv(path_output_low, sep=';', decimal='.', index_col='Unnamed: 0')
    df_output_high = read_csv(path_output_high, sep=';', decimal='.', index_col='Unnamed: 0')

    fig, ax = pyplot.subplots()
    for s, df, c in (('clairsemé', df_output_low, 'orange'), ('dense', df_output_high, 'red')):
        df.index = to_datetime(df.index)
        ax.bar(s, df['An'].sum() * CONV_CO2 / df['E'].sum(), facecolor=c)
    ax.set(ylabel='\n'.join(["Efficience de l'utilisation de l'eau", "$\mathregular{g_{CO_2}\/g^{-1}_{H_2O}}$"]))
    fig.tight_layout()

    pass


if __name__ == '__main__':
    path_data = Path(__file__).parent / 'data'
    # generate_grapevine_mtgs(path_dir_project=path_data)
    # %gui qt5

    # for is_low_area in (True, False):
    #     vine_mtg, output = run_hydroshoot(
    #         path_dir_preprocessed_inputs=path_data / 'preprocessed_inputs',
    #         is_low_leaf_area=is_low_area,
    #         is_show=True,
    #         path_output=path_data / f'output_{"low" if is_low_area else "high"}/time_series.csv')

    display_whole_plant(
        path_output_low=path_data / 'output_low/time_series.csv',
        path_output_high=path_data / 'output_high/time_series.csv',
        path_weather=path_data / 'weather.csv')

    display_mtg_properties(
        path_output_low=path_data / 'output_low',
        path_output_high=path_data / 'output_high')

    display_water_use_efficiency(
        path_output_low=path_data / 'output_low/time_series.csv',
        path_output_high=path_data / 'output_high/time_series.csv')
