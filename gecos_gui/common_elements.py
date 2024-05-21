import PySimpleGUI as Sg

theme = "LightGrey3"
Sg.set_options(text_justification='left', font=("Noto Sans", 11))
Sg.ChangeLookAndFeel(theme)

row_buttons = Sg.Column([
        [Sg.Button('Import data from Python', key='-BUTTONIMPORTPYTHON-', disabled=False),
         Sg.Button('Export data from Python', key='-BUTTONEXPORTPYTHON-', disabled=True),
         Sg.Button('Check Inputs', key='-BUTTONCHECK-', disabled=True),
         Sg.Button('Run/Check GeCos', key='-BUTTONRUN-', disabled=True),
         Sg.Button('Create Python Script', key='-BUTTONCREATESCRIPT-', disabled=True),
         Sg.Button('Visualize Results', key='-BUTTONVISRESULTS-', disabled=True),
         Sg.Button('View Log', key='-BUTTONVIEWLOG-', disabled=True)]
   ], pad=((90, 10), (0, 10)))

status_label = Sg.Text(text="Status: No data loaded",
                       key='-STATUS_TEXT-', size=(260, 1),background_color="lightblue", text_color="red",
                       font=('Expansiva', 10, 'bold'), justification='center')
# suggest_label = Sg.Text(text="Advice: Import data from Json file, python script or fill the form.",
#                         key='-SUGGEST_TEXT-', size=(60, 1))
inforun_label = Sg.Text(text="Advice: Import data from Json file, python script or fill the form.",
                        key='-INFORUN_TEXT-', size=(260, 1), background_color="lightblue",
                        font=('Expansiva', 10, 'bold'), justification='center')

python_script = Sg.Text('', key='-HIDEPYTHONSCRIPT-', visible=False)