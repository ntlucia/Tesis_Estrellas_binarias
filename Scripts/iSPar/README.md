# DataRedu

## Descripción: 
iSPar es un script en python que llama algunas de las librerias del programa iSpec, pero está modificado para
normalizar espectros de estrellas en secuencia principal y estrellas gigantes principalmante. Este código tiene una mayor
precisión en el canal rojo, ya que tiene una relación señal a ruido mucho mayor, sin embargo tiene buenos resultados para
el canal azul. Con este script también se pueden calcular los parámetros estelares, ya sea dandole un input o haciendo iteraciones y
comparaciones con el espectro sintético de la estrella.


**DataRedu** puede convertir recursivamente archivos .fit a .dat; esto es, todos los archivos de una carpeta o 
de carpetas separadas (cada carpeta con el nombre de la estrella).

#### Vemos que tenemos dos carpetas, una de ellas aplica la normalización que viene directa de la universidad de 
Hamburgo sobre todos los datos del tigre, pero es posible no usar esta normalización y hacerla con el código de iSPar
que está en la otra carpeta.

## Requisitos
* Tener instalado iSpec: Si no lo tiene instalado lo puede hacer [aquí](https://www.blancocuaresma.com/s/iSpec/manual/introduction)
* Python 2.7
* S.O. compatible con Unix.

## Cómo instalar
Guardar en la carpeta donde está el o los archivos que se quieren normalizar o calcular los parámetros, los archivos
"iSPar_v3.9.py" y "iSPar_Config.py"

## Cómo usar
1 Mover el archivo "initial_estimate_grid_ranges.tsv" a $ispec_dir/input/minigrid

2 Mover el archivo "iSPar_Line_Regions_Filtered_Combined.dat" a $ispec_dir/input/regions

3 Mover el archivo "telluric_and_absorption.txt" s $ispec/input/regions/strong_lines

4 Editar el archivo "iSPar_Config.py" según la estrella que se está analizando

5 Ejecute en la terminal dentro de la carpeta con los archivos a normalizar o calcular parámetros
.\iSPar_v3.9.py


