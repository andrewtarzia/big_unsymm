import pathlib


def xtb_path():
    return pathlib.Path(
        '/home/atarzia/anaconda3/envs/big_unsymm/bin/xtb'
    )


def liga_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/liga/'
    )


def meta_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/meta/'
    )


def figu_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/figures/'
    )


def cage_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/cages/'
    )


def pd24_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/m24_structures/'
    )


def manu_subs_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/cages/manual_extraction/'
    )


def subs_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/subsystems/'
    )


def calc_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/calculations/'
    )


def crest_path():
    return pathlib.Path(
        '/data/atarzia/projects/big_unsymm/software/crest'
    )


def crest_conformer_settings(solvent=None):
    return {
        'conf_opt_level': 'crude',
        'final_opt_level': 'extreme',
        'charge': 0,
        'no_unpaired_e': 0,
        'max_runs': 1,
        'calc_hessian': False,
        'solvent': solvent,
        'nc': 4,
        'etemp': 300,
        'keepdir': False,
        'cross': True,
        'md_len': None,
        'ewin': 5,
        'speed_setting': None,
    }


def gulp_path():
    return pathlib.Path(
        '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'
    )
