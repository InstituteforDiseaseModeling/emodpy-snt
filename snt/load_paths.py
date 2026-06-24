import os


def load_box_paths(user_path=None, country_name='Example'):
    if not user_path:
        user_path = os.path.expanduser('~')

    if country_name == 'Example':
        home_path = user_path
        data_path = os.path.join(user_path, 'Documents', 'emodpy-snt', 'data', 'example_files')
        project_path = os.path.join(user_path, 'Documents', 'emodpy-snt', 'data', 'example_files')
    elif country_name == 'SierraLeone':
        home_path = os.path.join(user_path, 'IDM Dropbox', 'Malaria Team Folder')
        data_path = os.path.join(home_path, 'data')
        project_path = os.path.join(home_path, 'projects', 'SierraLeone_hbhi')
    elif country_name == 'Burundi':
        home_path = os.path.join(user_path, 'IDM Dropbox', 'Malaria Team Folder')
        data_path = os.path.join(home_path, 'data')
        project_path = os.path.join(home_path, 'projects', 'burundi_hbhi', 'snt_2023')
    elif country_name == 'Nigeria':
        # home_path = os.path.join(user_path, 'Dropbox (IDM)', 'NU_collaboration')
        # data_path = os.path.join(home_path, 'hbhi_nigeria', 'snt_2022')
        # project_path = os.path.join(home_path, 'hbhi_nigeria', 'snt_2022')
        home_path = os.path.join(user_path, 'Gates Foundation Dropbox', 'Malaria Team Folder', 'projects', 'snt')
        data_path = os.path.join(home_path, 'Nigeria', 'snt_2026')
        project_path = os.path.join(home_path, 'Nigeria', 'snt_2026')
    return data_path, project_path


def load_emodpy_snt_path():
    user_path = os.path.expanduser('~')
    script_path = os.path.join(user_path, 'Documents', 'emodpy-snt')
    return script_path