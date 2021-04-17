# master-thesis
Minimum-width grid graph drawing algorithms

### Install libcairo2-dev
``` 
curl -L https://www.cairographics.org/releases/cairo-1.14.6.tar.xz -o cairo.tar.xz
tar -xf cairo.tar.xz && cd cairo-1.14.6
./configure --prefix=/usr/local --disable-dependency-tracking
make install
```

- On MacOS install pkg-config as well
``` brew install pkg-config ```


Run ```pip3 install -r requirements``` to install requirements.

To run project run ```python3 graph_generation/__main.py```