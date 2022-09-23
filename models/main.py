from json import dump
from json import load as load_json
from pathlib import Path
from pickle import load as load_pickle

from hydroshoot import architecture, display, io, initialisation, model
from openalea.mtg import traversal, mtg
from openalea.plantgl.all import Scene
from openalea.plantgl.all import surface as surf


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


if __name__ == '__main__':
    path_data = Path(__file__).parent / 'data'
    # generate_grapevine_mtgs(path_dir_project=path_data)
    # %gui qt5

    for is_low_area in (True, False):
        vine_mtg, output = run_hydroshoot(
            path_dir_preprocessed_inputs=path_data / 'preprocessed_inputs',
            is_low_leaf_area=is_low_area,
            is_show=True,
            path_output=path_data / f'output_{"low" if is_low_area else "high"}/time_series.csv')
