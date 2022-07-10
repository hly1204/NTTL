# Number Theory Template Library

NTTL 是一个只有头文件的数论相关的基础算法库, 要使用它只需要包含 `./src` 中的头文件即可.

## 测试

在本地测试需要使用比较高版本的 GCC 或者 MSVC, 我使用 GCC 12.1.0 进行测试.

```sh
cmake -S . -B build -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
cd build && make
cd test/unit_test/ && ctest
```
