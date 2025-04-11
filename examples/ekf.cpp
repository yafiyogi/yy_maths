/*

  MIT License

  Copyright (c) 2025 Yafiyogi

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "boost/numeric/ublas/assignment.hpp"
#include "fmt/format.h"

#include "yy_matrix_fmt.hpp"
#include "yy_matrix_util.hpp"
#include "yy_ekf.hpp"

using namespace yafiyogi;

using ekf = yy_maths::ekf;
using value_type = ekf::value_type;
using size_type = ekf::size_type;
using matrix = ekf::matrix;
using vector = ekf::vector;

constexpr size_type ekf_n = 2;
constexpr size_type ekf_m = 3;

void predict_update(ekf & filter,
                    const matrix & h,
                    const vector & z)
{
  filter.predict();

  vector hx{ekf_m};
  hx <<= filter.X(0), filter.X(1), filter.X(1);

  filter.update(z, h, hx);
}

int main()
{
  matrix H{ekf_m, ekf_n};
  H <<= 1, 0,
        0, 1,
        0, 1;

  matrix H1{ekf_m, ekf_n};
  H1 <<= 1, 0,
         0, 1,
         0, 0;

  matrix H2{ekf_m, ekf_n};
  H2 <<= 1, 0,
         0, 0,
         0, 1;

  vector z{ekf_m};

  {
    ekf ekf_state{ekf_m, ekf_n};

    z <<= 49.49, 22, 23;
    predict_update(ekf_state, H, z);
    fmt::print("1) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.59, 21, 24;
    predict_update(ekf_state, H, z);
    fmt::print("2) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));
  }

  fmt::print("\n");

  {
    ekf ekf_state{ekf_m, ekf_n};

    z <<= 49.49, 22, 23;
    predict_update(ekf_state, H1, z);
    fmt::print("1a) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.49, 22, 0.0;
    predict_update(ekf_state, H1, z);
    fmt::print("1a) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.49, 22, 23;
    predict_update(ekf_state, H1, z);
    fmt::print("1a) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.49, 22, 23;
    predict_update(ekf_state, H2, z);
    fmt::print("1b) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.59, 21, 24;
    predict_update(ekf_state, H1, z);
    fmt::print("2a) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));

    z <<= 49.59, 0.0, 24;
    predict_update(ekf_state, H2, z);
    fmt::print("2b) h=[{}] t=[{}]\n", ekf_state.X(0), ekf_state.X(1));
  }

  return 0;
}

// Study Temp
/*
  {"battery":100,"humidity":49.49,"last_seen":"2025-01-07T11:19:15.194Z","linkquality":127,"temperature":22,"voltage":3000}
  {"battery":100,"humidity":49.59,"last_seen":"2025-01-07T11:21:27.505Z","linkquality":131,"temperature":22,"voltage":3000}
  {"battery":100,"humidity":49.59,"last_seen":"2025-01-07T11:23:26.631Z","linkquality":127,"temperature":22,"voltage":3000}
  {"battery":100,"humidity":49.59,"last_seen":"2025-01-07T11:23:26.651Z","linkquality":131,"temperature":22,"voltage":3000}
  {"battery":100,"humidity":49.59,"last_seen":"2025-01-07T11:24:15.294Z","linkquality":131,"temperature":22.05,"voltage":3000}
  {"battery":100,"humidity":49.57,"last_seen":"2025-01-07T11:26:27.601Z","linkquality":134,"temperature":22.05,"voltage":3000}
  {"battery":100,"humidity":49.57,"last_seen":"2025-01-07T11:29:15.416Z","linkquality":127,"temperature":22.09,"voltage":3000}
  {"battery":100,"humidity":49.76,"last_seen":"2025-01-07T11:31:27.727Z","linkquality":123,"temperature":22.09,"voltage":3000}
  {"battery":100,"humidity":49.76,"last_seen":"2025-01-07T11:34:15.540Z","linkquality":120,"temperature":22.12,"voltage":3000}
*/

// AirQM
/*
  {"humidity":44,"last_seen":"2025-01-07T11:19:56.826Z","linkquality":131,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":99}
  {"humidity":44,"last_seen":"2025-01-07T11:20:56.485Z","linkquality":131,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:21:48.197Z","linkquality":131,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:23:47.913Z","linkquality":127,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":103}
  {"humidity":44,"last_seen":"2025-01-07T11:24:36.472Z","linkquality":127,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:24:55.916Z","linkquality":127,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:25:47.630Z","linkquality":123,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:25:55.975Z","linkquality":123,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
  {"humidity":44,"last_seen":"2025-01-07T11:26:36.189Z","linkquality":123,"pm25":1,"temperature":23,"update":{"installed_version":16777233,"latest_version":16777233,"state":"idle"},"update_available":null,"voc_index":102}
*/
