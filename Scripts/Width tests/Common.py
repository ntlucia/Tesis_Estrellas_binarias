import yaml

__config = None
__Spectra = None

a=1
def config():
    global __config
    if not __config: #Solo accede a disco si no ha sido cargado el archivo de config.yaml
        with open("config.yaml",mode='r') as f:
            __config = yaml.load(f)
    return __config

def Spectra():
    global __Spectra
    if not __Spectra: 
        with open("Spectra_routes.txt",mode='r') as f:
            __Spectra = f.readlines()
    return __Spectra[0]
