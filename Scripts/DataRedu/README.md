# DataRedu

## Descripción: 
DataRedu es un script en bash que llama a una versión modificada por mí del programa en python 
con_to_ascii.py diseñado por el grupo de investigación de atmósferas estelares del convenio internacional 
entre la Universidad de Hamburgo (Alemania), la Universidad de Guanajuato (México) y la Universidad de Liegè 
(Bélgica); para convertir los archicos .fit de los espectros tomados por TIGRE-HEROS a archivos .dat


Con la versión modificada de con_to_ascii.py, todos los .fits del directorio serán convertidos recursivamente
 a .dat para ser leídos correctamente en iSpec. Los .dat salientes tendrán la siguiente estructura de columnas 
de izquierda a derecha: Wavelength (nm), Flux (erg/cm^2/s) y Flux Error (erg/cm^2/s)

**DataRedu** puede convertir recursivamente archivos .fit a .dat; esto es, todos los archivos de una carpeta o 
de carpetas separadas (cada carpeta con el nombre de la estrella).

#### Vemos que tenemos dos carpetas, una de ellas aplica la normalización que viene directa de la universidad de 
Hamburgo sobre todos los datos del tigre, pero es posible no usar esta normalización y hacerla con el código de iSPar
que está en la otra carpeta.

## Requisitos

* Python 2.7
* S.O. compatible con Unix.

## Cómo instalar

Desde la carpeta donde vayas a guardar el script, ejecutar en la terminal

./DataRedu_Compiler.sh

Esto generará 2 archivos: DataRedu.sh y DataRedu_iSpec.sh


## Cómo usar
1 El nombre de la carpeta que contiene los .fit debe ser el mismo que el nombre de la estrella 
en el archivo, por ejemplo:

Si el archivo .fit tiene de nombre: zetaAur_eclipse_B_2019_12_18_19_27_58.fits

La carpeta contenedora debe llamarse: zetaAur_eclipse

La anterior condición es necesaria para el correcto funcionamiento del script.

2 Desde la carpeta general donde se tiene las distintas subcarpetas de cada estrella, copiar y 
ejecutar DataRedu.sh, esto copiará DataRedu_iSpec.sh en cada subcarpeta y ejecutará la reducción para
 todas las estrellas pasando de .fit a .dat

El script también obtendrá información relevante del canal R de HEROS, como resolución media (Resolution), 
relación señal a ruido (SNR) y el valor de escala de la función de respuesta instrumental (IRF), este último 
es necesario para igualar la respuesta del canal R y el canal B de HEROS. La anterior información se guarda 
en el archivo <star_name>_INFO.dat

Nota: Si se desea repetir la reducción únicamente dentro de una carpeta de estrella en particular, solo basta con 
ejecutar DataRedu_iSpec.sh.