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
sudo apt install libgmp3-dev
sudo apt install libc6
sudo apt upgrade libstdc++6
```

Install Clang 15:

```console
sudo apt-key adv --keyserver keyserver.ubuntu.com:80 --recv-keys A4B469963BF863CC
sudo apt install lsb-release wget software-properties-common
sudo bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
```

Compile program:

```console
sudo clang++-15 -std=c++20 -O3 -m64 congruent_numbers.cpp -o congruent_numbers -lpthread -lstdc++
```
or
```console
sudo /usr/bin/g++-11 -std=c++20 -O3 -m64 congruent_numbers.cpp -o congruent_numbers -lpthread -lstdc++
```

Run program:

```console
sudo ./congruent_numbers
```

Obtain a status:

```console
tail -c 2000000 ec_pq_find_rational_xy.log | grep "Total PQs"
```

The command above yields for example:

```console
--------------- Total PQs 40065 (last 5,45263), Solved 6963 (17.4%), Worst 3,11 (1 blocks) ---------------
--------------- Total PQs 40066 (last 227,997), Solved 6963 (17.4%), Worst 3,11 (1 blocks) ---------------
--------------- Total PQs 40067 (last 17,13313), Solved 6963 (17.4%), Worst 3,11 (1 blocks) ---------------
--------------- Total PQs 40068 (last 89,2543), Solved 6963 (17.4%), Worst 3,11 (1 blocks) ---------------
```

***Note:*** Compile and run the program `congruent_numbers_ab.cpp` analoguosly to `congruent_numbers.cpp`.

## Compiling and running under Windows
Start MSYS2 MinGW Clang x64

Compile program:

```console
g++ -std=c++20 -O3 -m64 /c/Users/esultano/git/elliptic_curves/cpp/congruent_numbers.cpp -o /c/Users/esultano/git/elliptic_curves/cpp/congruent_numbers -lpthread -lstdc++ -lm
```

Run program:

```console
/c/Users/esultano/git/elliptic_curves/cpp/congruent_numbers.exe
```
