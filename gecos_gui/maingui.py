import PySimpleGUI as Sg
import os
from gecos_gui.main_tabgroup import tabgroup_layout_col
from gecos_gui.menu_layout import menu_layout
from gecos_gui.main_figure_layout import image_layout, lizard_gif
from gecos_gui.events import waiting_for_events
from gecos_gui.common_elements import row_buttons, theme, status_label, python_script, inforun_label


def main_gui_app():

    # MAIN WINDOW LOOP
    Sg.ChangeLookAndFeel(theme)
    layout = [[menu_layout], [image_layout, tabgroup_layout_col],
              [status_label], [inforun_label], [python_script]]

    window = Sg.Window('GeCos GUI', layout, size=(1650, 890), location=(100, 60),
                       finalize=True, return_keyboard_events=True, force_toplevel=False)

    window['-DOCK_URL-'].bind('<Enter>', '+MOUSE OVER+')
    window['-DOCK_URL-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-CONF_URL-'].bind('<Enter>', '+MOUSE OVER+')
    window['-CONF_URL-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-LINK_AUTHOR-'].bind('<Enter>', '+MOUSE OVER+')
    window['-LINK_AUTHOR-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-LINK_DOCS-'].bind('<Enter>', '+MOUSE OVER+')
    window['-LINK_DOCS-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-LINK_GITHUB-'].bind('<Enter>', '+MOUSE OVER+')
    window['-LINK_GITHUB-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-GAUSSIAN16_URL-'].bind('<Enter>', '+MOUSE OVER+')
    window['-GAUSSIAN16_URL-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-GAUSSIAN16_EXTRAINFO_URL-'].bind('<Enter>', '+MOUSE OVER+')
    window['-GAUSSIAN16_EXTRAINFO_URL-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-GAUSSIAN16_EXTRAINFO-'].bind('<Enter>', '+MOUSE OVER+')
    window['-GAUSSIAN16_EXTRAINFO-'].bind('<Leave>', '+MOUSE AWAY+')
    window['-INPUT_OTHERSOLV_PROP-'].bind('<Enter>', '+MOUSE OVER+')
    window['-INPUT_OTHERSOLV_PROP-'].bind('<Leave>', '+MOUSE AWAY+')

    while True:

        # event, values = window.read(timeout=125)
        event, values = window.read()
        if event == Sg.WIN_CLOSED or event == "Exit":
            script_path = window['-HIDEINPUTSCRIPT-'].get()
            directory = os.path.dirname(script_path)+"/"
            if os.path.exists(directory+"tmp.sh"):
                os.remove(directory+"tmp.sh")
            break

        window.Element('-IMAGE_LIZARD-').UpdateAnimation(lizard_gif,  time_between_frames=50)

        waiting_for_events(window, event, values)

    window.close()


if __name__ == "__main__":
    main_gui_app()
