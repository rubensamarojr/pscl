# Modified project mixing libigl, fcpw and fpdc

    git clone --recursive https://github.com/rubensamarojr/fcpw-libigl-example/

## Compile

Go to the folder **fcpw-libigl-example**

    cd fcpw-libigl-example

Run the command cmake

    cmake .

Execute make

    make

This should create the executable `fcpw-libigl-fpdc-example.exe`

## Run

### Windows
    fcpw-libigl-fpdc-example.exe [path-to-mesh] 

    fcpw-libigl-fpdc-example.exe [path-to-mesh] [path-to-queries].dmat

### Linux
    ./fcpw-libigl-fpdc-example [path-to-mesh] 

    ./fcpw-libigl-fpdc-example [path-to-mesh] [path-to-queries].dmat

Will output 

    - `Q.dmat`  #Q by 3 list of query points
    - `I.dmat`  #I list of indices into triangles
    - `UV.dmat`  #Q by 2 list of barycentric coordinates

