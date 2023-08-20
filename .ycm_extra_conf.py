import os

# Base flags
flags = [
    '-x',
    'c++',
    '-Wall',
    '-Wextra',
    '-Werror',
]


# Exclude directories with specific name
EXCLUDE_DIRS = [
    "vs2008",
]


def MoveDirectoryUpUntil(directory, list_required=[]) -> str:
    required = list_required.copy()
    current_content = os.listdir(directory)

    for current_object in current_content:
        if current_object in required: required.remove(current_object)

    if required != []:
        directory_up = os.path.split(directory)[0]
        if (len(directory_up) <= 6) | (".git" in current_content): # i.e. /home/ or highest up in project
            print(f"WARNING: No directory found containing: {list_required}")
            return None
        directory = MoveDirectoryUpUntil(directory_up, list_required)
    return directory


def GetAllSubDirectories(directory):
    ret = []
    for dirpath, dirnames, filenames in os.walk(directory):
        dirnames.append('')
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
        for dirname in dirnames:
            ret.append(os.path.join(dirpath, dirname))
        if '' in dirnames: dirnames.remove('')
    return ret


def GetProjectDir(directory):
    project_dir = MoveDirectoryUpUntil(directory, ['.git'])
    return project_dir


def AddAllDirectories(flags, working_directory):
    proj_dir = GetProjectDir(working_directory)

    folders = GetAllSubDirectories(proj_dir)
    if folders != []:
        for fold in folders:
            flags.append(f'-I{fold}')

    return flags


# Get the relevant directories and add the flags
working_directory = os.path.abspath('.')
flags = AddAllDirectories(flags, working_directory)

def Settings( **kwargs ):
    return {'flags': flags}
