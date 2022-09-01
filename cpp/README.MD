## Arty's programm to search integer solutions of an elliptic curve
With the greatful help of [Arty](https://chat.stackoverflow.com/users/941531/arty) we could demonstrate how `C++` beats the performance and precision of other approaches involving `Python` along with `JIT`, `CUDA`.

On Linux do the following:

#### compile

```console
g++-11 -std=c++20 -O3 -m64 prog.cpp -o prog -lpthread -lgmp
```

#### run

```console
./prog 0 100000000
```

## Compiling and running under Windows
On Windows do the following:

#### install libgmp

```console
pacman -S mingw-w64-x86_64-gmp
```

#### compile

```console
g++ -std=c++20 -O3 -m64 /c/Users/esultano/git/elliptic_curves/cpp/prog.cpp -o /c/Users/esultano/git/elliptic_curves/cpp/prog -lpthread -lgmp
```

#### run

```console
/c/Users/esultano/git/elliptic_curves/cpp/prog 0 100000000
```

## Arty's programm for congruent numbers

Install prerequisites:

```console
apt install sudo
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt upgrade
sudo apt install gcc
sudo apt install g++-11
sudo apt install libc6
sudo apt upgrade libstdc++6
```

Compile program:

```console
g++-11 -std=c++20 -O3 -m64 congruent_numbers.cpp -o congruent_numbers -lpthread
```

Run program:

```console
./congruent_numbers
```
