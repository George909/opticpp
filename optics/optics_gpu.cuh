#include <thrust/complex.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>

#include <complex>
#include <vector>
#include <array>

namespace optics_gpu {

    using complex = std::complex<double>;    
    using vec3_complex = std::array<complex, 3>;
    using vec_vec3_complex = std::vector<vec3_complex>;
    using vec_complex = std::vector<complex>;

    using dev_complex = thrust::complex<double>;    
    using dev_vec_complex = thrust::device_vector<dev_complex>;

    struct dev_vec3 {
        dev_complex x;
        dev_complex y;
        dev_complex z;
    };

    using dev_vec_vec3 = thrust::device_vector<dev_vec3>;

    struct rw_info {
        double na;
        double lambda;
        double out_w;
        double out_h;
        double in_w;
        double in_h;
    };

    struct linear {
        __host__ __device__
        dev_vec3 operator()(const double& phi, const double& theta) const {
            double cos_theta = cos(theta);
            double cos_phi = cos(phi);
            return {1 + cos_phi * cos_phi * (cos_theta - 1),
                    sin(phi) * cos_phi * (cos_theta - 1),
                    -cos_phi * sin(theta)};
        };
    };

    struct rcircular {
        __host__ __device__
        dev_vec3 operator()(const double& phi, const double& theta) const {
            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            double cos_theta = cos(theta);
            complex I = complex(0.,1.);

            return {M_SQRT1_2 * (1 + cos_phi * cos_phi * (cos_theta - 1) + I * (sin_phi * cos_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (sin_phi * cos_phi * (cos_theta - 1) + I * (1 + sin_phi * sin_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (-sin(theta) * (cos_phi + I * sin_phi))};
        };
    };

    struct lcircular {
        __host__ __device__
        dev_vec3 operator()(const double& phi, const double& theta) const {
            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            double cos_theta = cos(theta);
            complex I = complex(0.,1.);

            return {M_SQRT1_2 * (1 + cos_phi * cos_phi * (cos_theta - 1) - I * (sin_phi * cos_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (sin_phi * cos_phi * (cos_theta - 1) - I * (1 + sin_phi * sin_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (-sin(theta) * (cos_phi - I * sin_phi))};
        };
    };

    struct radial {
        __host__ __device__
        dev_vec3 operator()(const double& phi, const double& theta) const {
            double cos_theta = cos(theta);

            return {cos(phi) * cos_theta,
                    sin(phi) * cos_theta,
                    -sin(theta)};
        }
    };

    struct azimutal {
        __host__ __device__
        dev_vec3 operator()(const double& phi, const double& theta) const {
            return {sin(phi),
                    -cos(theta),
                    0.};
        }
    };

    template<typename Functor>
    void rw(vec_complex& input, size_t w, size_t h, 
            rw_info& inf,vec_vec3_complex& output,
            Functor func = {}) {
        dev_vec_complex dev_input = input;
        dev_vec_vec3 dev_output;

        dev_complex I = dev_complex(0.,1.);
        double out_w_ratio = inf.out_w / w ;
        double out_h_ratio = inf.out_h / h;
        double in_w_ratio = inf.in_w / w;
        double in_h_ratio = inf.in_h / h;
        double k = 2 * M_PI / inf.lambda;
        double f = inf.in_w / (2. * tan(inf.na));
        double lambda = inf.lambda;

        auto sum_up ( __device__ [=, &dev_input](size_t out_index){
            double x_out = (out_index % w + 1 - w / 2.) * out_w_ratio;
            double y_out = (out_index / w + 1 - h / 2.) * out_h_ratio;
            double phi_out = atan2(y_out, x_out);
            double r_out = sqrt(x_out * x_out + y_out * y_out);
            return __device__ [=, &dev_input](size_t in_index, dev_complex c) -> dev_vec3 {
                double x_in = (in_index % w + 1 - w / 2.) * in_w_ratio;
                double y_in = (in_index / w + 1 - h / 2.) * in_h_ratio;
                double r_in = sqrt(x_in * x_in + y_in * y_in);
                double phi_in = atan2(y_in, x_in);
                double theta = atan2(r_in, f);
                double sin_theta = sin(theta);
                dev_vec3 p = func(phi_in, theta);
                dev_complex coef = -sqrt(cos(theta)) * c
                               * exp(I * k * r_out * sin_theta * cos(phi_in - phi_out))
                               * sin_theta * I * f / lambda;

                return {coef * p.x, coef * p.y, coef * p.z};
            };
        } );
        
        auto plus_vec3 (__host__ __device__ 
                        [](const dev_vec3& r, const dev_vec3& l) -> dev_vec3 {
                            return {r.x + l.x, r.y + l.y, r.z + l.z};
        } );

        auto rw_tr (__host__ __device__ [=,&dev_input](size_t out_index){
            return thrust::transform_reduce(thrust::counting_iterator<size_t>{0},
                                            thrust::counting_iterator<size_t>{dev_input.size()},
                                            dev_input.begin(),
                                            sum_up(out_index),
                                            dev_vec3{},
                                            plus_vec3);
        });
        
        thrust::transform(thrust::counting_iterator<size_t>{0},
                          thrust::counting_iterator<size_t>{input.size()},
                          dev_output.begin(),
                          rw_tr);
        
        //output = std::reinterpret_cast<vec_vec3_complex>(dev_output);
    }
}
