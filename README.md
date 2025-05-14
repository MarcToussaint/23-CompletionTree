# 23-CompletionTree

Source code for this paper:
https://www.user.tu-berlin.de/mtoussai/24-CompletionTrees/

Try `make` src/

Dependencies need to be installed, e.g. as for the https://github.com/MarcToussaint/robotic package. I think this might be sufficient:

```
sudo apt install --yes \
  g++ clang make gnupg cmake git wget libstdc++-14-dev \
  liblapack-dev libf2c2-dev libqhull-dev libeigen3-dev \
  libjsoncpp-dev libyaml-cpp-dev libhdf5-dev libpoco-dev libboost-system-dev portaudio19-dev libusb-1.0-0-dev \
  libx11-dev libglu1-mesa-dev libglfw3-dev libglew-dev freeglut3-dev libpng-dev libassimp-dev
  
wget https://github.com/MarcToussaint/rai/raw/refs/heads/marc/_make/install.sh; chmod a+x install.sh
./install.sh libccd
./install.sh fcl
./install.sh libann
```
