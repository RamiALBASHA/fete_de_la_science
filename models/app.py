# Adapted from https://shinylive.io/py/examples/#orbit-simulation

from pathlib import Path

from matplotlib import pyplot
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_widget

import sim_funcs

pyplot.ioff()

path_root = Path(__file__).parent
path_data = path_root / 'data'
path_static = path_root / "static"


def panel_box(*args, **kwargs):
    return ui.div(
        ui.div(*args, class_="card-body"),
        **kwargs,
        class_="card mb-3",
    )


app_ui = ui.page_fluid(
    {"class": "p-4"},
    ui.row(
        ui.column(
            2,
            ui.img(src="logos_mesr_occitanie.png", style="width: 100%; max-width: 250px;"),
            ui.img(src="logo_fds_noir_rouge.png",
                   style="width: 50%; max-width: 250px; vertical-align:center, horizontal-align:center"),

        ),
        ui.column(8, ui.h1("Quelle conduite de la vigne face au changement climatique ?",
                           style="font-size: 30px; color: Navy; text-align: center; font-weight: bold;"
                                 "font-family:Arial Narrow; background-color:  #f8f9f9 ")),
        ui.column(
            2,
            ui.img(src="logo_lepse.jpg", style="width: 100%; max-width: 250px;"),
            ui.img(src="logo_instant_science.png", style="width: 50%; max-width: 250px;")
        ),
    ),
    ui.row(
        ui.column(
            5,
            panel_box(
                ui.row(
                    ui.column(6,
                              ui.input_select(
                                  "canopy_1",
                                  "Conduite 1",
                                  {
                                      "vsp": "Espalier bas",
                                      "gdc": "Rideau simple",
                                      "lyre": "Lyre",
                                      "cordon": "Cordon libre",
                                  },
                              ),
                              ),
                    ui.column(6,
                              ui.input_select(
                                  "canopy_2",
                                  "Conduite 2",
                                  {
                                      "vsp": "Espalier bas",
                                      "gdc": "Rideau simple",
                                      "lyre": "Lyre",
                                      "cordon": "Cordon libre",
                                  },
                              ),
                              ),
                ),
                ui.row(
                    ui.column(6,
                              # ui.output_ui(id="display_canopy1"),
                              output_widget(id="display_3d_canopy1")
                              ),
                    ui.column(6,
                              # ui.output_ui(id="display_canopy2"),
                              output_widget(id="display_3d_canopy2")
                              )
                ),
                ui.input_action_button("compare", "Comparer", class_="btn-primary w-100"),
            ),
            panel_box(
                ui.input_numeric(id="hour", label="Heure simulée", value=12, min=9, max=15, step=1),
                ui.input_select(
                    id="mtg_prop",
                    label="Propriété à afficher",
                    choices={
                        "gs": "Conductance stomatique",
                        "An": "Photosynthèse",
                        "psi_head": "Potentiel Xylémien",
                        "Tlc": "Température",
                        "Na": "Teneur en azote",
                        "E": "Transpiration",
                        "gb": "Conductance Couche lim.",
                    }
                ),
                ui.input_action_button(id="plot_mtg", label="Comparer", class_="btn-primary w-100"),
            ),
            panel_box(
                ui.h2("Simulation à l'échelle de la feuille",
                      style="font-size: 20px; color: Navy; text-align: center"),
                ui.output_plot("display_mtg_property", height='200px')
            )
        ),
        ui.column(
            4,
            panel_box(
                ui.h2("Simulation à l'échelle de la canopée",
                      style="font-size: 20px; color: Navy; text-align: center"),
                ui.output_plot("display_whole_plant", height='500px'),
                ui.output_plot("display_water_use_efficiency", height='300px'),
            ),
        ),
        ui.column(
            3,
            panel_box(
                ui.row(
                    ui.column(8,
                              ui.input_select(
                                  id="canopy_user",
                                  label="Conduite 3",
                                  choices={
                                      "none": "aucun",
                                      "vsp": "Espalier bas",
                                      "gdc": "Rideau simple",
                                      "lyre": "Lyre",
                                      "cordon": "Cordon libre",
                                  },
                                  selected=None,
                              ),
                              ),
                    ui.column(4,
                              ui.output_ui("display_canopy_user"),
                              )
                ),
                ui.input_slider("gs_min", "Conductance stomatique minimale\n(mmol/m2/s)", 0, 40, value=20),
                ui.input_slider("psi_gs50", "Pression à laquelle les stomates ferment à 50% (MPa)", -4, 0, value=-1,
                                step=0.1),
                ui.input_slider("n", "Sensibilité au déficit hydrique", 0, 5, value=4, step=0.1),
                ui.input_action_button("run", "Lancer", class_="btn-primary w-100"),
                ui.output_plot("display_stomatal_response"),
            ),
        )
    ),
)


def server(ui_input, ui_output, session):
    @ui_output
    @render.ui
    def display_canopy1():
        return ui.tags.img(src=f'icone_{ui_input.canopy_1()}.png', height='150px')

    @ui_output
    @render.ui
    def display_canopy2():
        return ui.tags.img(src=f'icone_{ui_input.canopy_2()}.png', height='150px')

    @ui_output
    @render.ui
    def display_canopy_user():
        return ui.tags.img(src=f'icone_{ui_input.canopy_user()}.png', height='100px')

    @ui_output
    @render.plot(alt="Simulation à l'échelle de la plante")
    @reactive.event(ui_input.compare, ui_input.run, ignore_init=True, ignore_none=False)
    def display_whole_plant():
        return sim_funcs.plot_whole_plant_gas_exchange(
            path_ref=path_data / 'output_ref',
            path_user=path_data / 'output_user',
            training_system_1=ui_input.canopy_1(),
            training_system_2=ui_input.canopy_2(),
            training_system_user=ui_input.canopy_user(),
            path_weather=path_data / 'weather.csv')

    @ui_output
    @render.plot
    @reactive.event(ui_input.compare, ui_input.run, ignore_init=True, ignore_none=False)
    def display_water_use_efficiency():
        return sim_funcs.plot_water_use_efficiency(
            path_ref=path_data / 'output_ref',
            path_user=path_data / 'output_user',
            training_system_1=ui_input.canopy_1(),
            training_system_2=ui_input.canopy_2(),
            training_system_user=ui_input.canopy_user())

    pass

    @ui_output
    @render.plot
    def display_stomatal_response():
        return sim_funcs.plot_stomatal_reduction_coef(
            psi50=ui_input.psi_gs50(),
            gs0=ui_input.gs_min() / 1000,
            n=ui_input.n())

    @ui_output
    @render.plot
    @reactive.event(ui_input.plot_mtg, ignore_init=True, ignore_none=False)
    def display_mtg_property():
        return sim_funcs.plot_mtg_property(
            path_ref=path_data / 'output_ref',
            path_user=path_data / 'output_user',
            training_system_1=ui_input.canopy_1(),
            training_system_2=ui_input.canopy_2(),
            training_system_user=ui_input.canopy_user(),
            hour=ui_input.hour(),
            mtg_property=ui_input.mtg_prop()
        )

    @reactive.Effect
    @reactive.event(ui_input.run, ignore_init=True, ignore_none=False)
    def update_model_params():
        sim_funcs.update_params(
            path_default_params=path_data / 'params_default.json',
            path_user_params=path_data / 'params.json',
            stomatal_conductance_min=ui_input.gs_min(),
            stomatal_psi50=ui_input.psi_gs50(),
            steepness=ui_input.n())

    @reactive.Effect
    @reactive.event(ui_input.run, ignore_init=True, ignore_none=False)
    def run_hydroshoot():
        path_output = path_data / 'output_user'
        path_output.mkdir(parents=True, exist_ok=True)
        sim_funcs.run_hydroshoot(
            path_project=path_data,
            path_dir_preprocessed_inputs=path_data / 'preprocessed_inputs',
            training_sys=ui_input.canopy_user(),
            is_show=False,
            path_output=path_output / 'time_series.csv')

    @ui_output
    @render_widget
    def display_3d_canopy1():
        return sim_funcs.show_3d_canopy(canopy_name=ui_input.canopy_1(), is_complete_canopy=True)

    @ui_output
    @render_widget
    def display_3d_canopy2():
        return sim_funcs.show_3d_canopy(canopy_name=ui_input.canopy_2(), is_complete_canopy=True)


app = App(ui=app_ui, server=server, debug=True, static_assets=path_static)
