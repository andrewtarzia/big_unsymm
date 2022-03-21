from pymongo import MongoClient
from getpass import getpass


def xtb_path():
    return '/home/atarzia/anaconda3/envs/big_unsymm/bin/xtb'


def liga_path():
    return '/data/atarzia/projects/big_unsymm/liga/'


def meta_path():
    return '/data/atarzia/projects/big_unsymm/meta/'


def cage_path():
    return '/data/atarzia/projects/big_unsymm/cages/'


def crest_path():
    raise NotImplementedError()
    return '/home/atarzia/software/crest/crest'


def gulp_path():
    raise NotImplementedError()
    return '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'


def mongo_client():
    raise NotImplementedError()
    user = 'atarzia'
    password = getpass()
    return MongoClient(
        f'mongodb+srv://{user}:{password}@cluster0.32as5.mongodb.net'
        '/lc?retryWrites=true&w=majority'
    )
