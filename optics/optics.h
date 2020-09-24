#ifndef OPTICS_H
#define OPTICS_H

#include <concepts>
#include <vector>
#include <array>
#include <complex>
#include <valarray>
#include <algorithm>
#include <execution>

#include <boost/iterator/counting_iterator.hpp>

namespace optics {
    using namespace std::complex_literals;

    using complex = std::complex<double>;
    using vec_complex = std::vector<complex>;
    using vec3_complex = std::array<complex, 3>;
    using vec_vec3_complex = std::vector<vec3_complex>;

    struct gauss_ampl {
        double sigma;
        double __sigma2;
        gauss_ampl() = default;
        gauss_ampl(double sigma) 
            : sigma{sigma}
            , __sigma2{sigma * sigma} {};
        double operator() (const double x, const double y) const {
            return std::exp(-(x*x + y*y) / __sigma2);
        };
    };

    struct some_ampl {
        double sin_sigma;
        double sin_alpha;
        some_ampl() = default;
        some_ampl(double sin_sigma, double sin_alpha) 
            : sin_sigma{sin_sigma}, sin_alpha{sin_alpha} {};
        double operator()(const double x, const double y) const {
            double angle =  std::abs(std::atan2(y, x));
            double sin_angle = std::sin(angle);
            return std::exp(-sin_angle * sin_angle / sin_sigma / sin_sigma) * sin_angle / sin_alpha;
        };
    };

    struct vortex_phase {
        double m, n;
        vortex_phase() : m{}, n{1} {};
        vortex_phase(double m) : m{m}, n{1} {};
        vortex_phase(double m, double n) : m{m}, n{n} {};
        double operator() (const double x, const double y) const {
            double angle = std::atan2(y, x) + M_PI;
            return m * std::pow(angle, n);
        }
    };

    struct some_phase {
        double a, b, n, m;
        some_phase() = default;
        some_phase(double a, double b, double n, double m)
            : a{a}, b{b}, n{n}, m{m} {};
        double operator() (const double x, const double y) const {
            return a * std::pow(x, n) + b * std::pow(y, m);
        }; 
    };

    template<typename Ampl = gauss_ampl, 
             typename Phase = vortex_phase>
    vec_complex beam(size_t w, size_t h, Ampl ampl = {}, Phase phase = {}) {
        vec_complex field(w*h);
        std::transform(std::execution::par,
                       boost::counting_iterator<size_t>{0},
                       boost::counting_iterator<size_t>{w*h},
                       field.begin(),
                       [=](const size_t& index) -> complex {
                           double x = index % w + 1 - w / 2.;
                           double y = index / w + 1 - h / 2.;
                           return ampl(x,y) * std::exp(1i * phase(x,y));
                       });

        return field;
    }

    struct rw_info {
        double na;
        double lambda;
        double out_w;
        double out_h;
        double in_w;
        double in_h;
    };

    struct linear {
        vec3_complex operator()(const double& phi, const double& theta) const {
            double cos_theta = std::cos(theta);
            double cos_phi = std::cos(phi);
            return {1 + cos_phi * cos_phi * (cos_theta - 1),
                    std::sin(phi) * cos_phi * (cos_theta - 1),
                    -cos_phi * std::sin(theta)};
        };
    };

    struct rcircular {
        vec3_complex operator()(const double& phi, const double& theta) const {
            double sin_phi = std::sin(phi);
            double cos_phi = std::cos(phi);
            double cos_theta = std::cos(theta);

            return {M_SQRT1_2 * (1 + cos_phi * cos_phi * (cos_theta - 1) + 1i * (sin_phi * cos_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (sin_phi * cos_phi * (cos_theta - 1) + 1i * (1 + sin_phi * sin_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (-std::sin(theta) * (cos_phi + 1i * sin_phi))};
        };
    };

    struct lcircular {
        vec3_complex operator()(const double& phi, const double& theta) const {
            double sin_phi = std::sin(phi);
            double cos_phi = std::cos(phi);
            double cos_theta = std::cos(theta);

            return {M_SQRT1_2 * (1 + cos_phi * cos_phi * (cos_theta - 1) - 1i * (sin_phi * cos_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (sin_phi * cos_phi * (cos_theta - 1) - 1i * (1 + sin_phi * sin_phi * (cos_theta - 1))),
                    M_SQRT1_2 * (-std::sin(theta) * (cos_phi - 1i * sin_phi))};
        };
    };

    struct radial {
        vec3_complex operator()(const double& phi, const double& theta) const {
            double cos_theta = std::cos(theta);

            return {std::cos(phi) * cos_theta,
                    std::sin(phi) * cos_theta,
                    -std::sin(theta)};
        }
    };

    struct azimutal {
        vec3_complex operator()(const double& phi, const double& theta) const {
            return {std::sin(phi),
                    -std::cos(theta),
                    0.};
        }
    };

    template<typename T>
    concept isPolarization = requires (T& t,const double& a,const double& b) {
        {t(a,b)} -> std::same_as<vec3_complex>;
    };

    template<isPolarization Functor>
    void richards_wolf(const vec_complex& input, rw_info& inf,
                       size_t h, size_t w, vec_vec3_complex& output,
                       Functor functor = {}){
        double out_w_ratio = inf.out_w / w ;
        double out_h_ratio = inf.out_h / h;
        double in_w_ratio = inf.in_w / w;
        double in_h_ratio = inf.in_h / h;
        double k = 2 * M_PI / inf.lambda;
        double f = std::min(inf.in_h, inf.in_w) / (2. * std::tan(inf.na));
        double lambda = inf.lambda;

        auto sum_up ([=, &input](const size_t out_index){
            double x_out = (out_index % w + 1 - w / 2.) * out_w_ratio;
            double y_out = (out_index / w + 1 - h / 2.) * out_h_ratio;
            double phi_out = std::atan2(y_out, x_out);
            double r_out = std::sqrt(x_out * x_out + y_out * y_out);
            return [=, &input](const size_t in_index, const complex c) -> vec3_complex {
                double x_in = (in_index % w + 1 - w / 2.) * in_w_ratio;
                double y_in = (in_index / w + 1 - h / 2.) * in_h_ratio;
                double r_in = std::sqrt(x_in * x_in + y_in * y_in);
                double phi_in = std::atan2(y_in, x_in);
                double theta = std::atan2(r_in, f);
                double sin_theta = std::sin(theta);
                vec3_complex p = functor(phi_in, theta);
                complex coef = std::sqrt(std::cos(theta)) * c
                               * std::exp(1i * k * r_out * sin_theta * std::cos(phi_in - phi_out))
                               * sin_theta;

                return {coef * p[0], coef * p[1], coef * p[2]};
            };
        });

        auto plus_vec3_complex ([](const vec3_complex& r, const vec3_complex& l) -> vec3_complex {
            return {r[0] + l[0], r[1] + l[1], r[2] + l[2]};
        });

        auto rw_tr ( [=,&input](const size_t out_index){
            auto vec = std::transform_reduce(boost::counting_iterator<size_t>{0},
                                             boost::counting_iterator<size_t>{input.size()},
                                             input.begin(),
                                             vec3_complex{0,0,0},
                                             plus_vec3_complex,
                                             sum_up(out_index));
            std::transform(vec.begin(), vec.end(), vec.begin(),[=](const auto& c) {
                                                                     return -1i * c * f / lambda;});
            return vec;
        } );

        std::transform(std::execution::par,
                       boost::counting_iterator<size_t>{0},
                       boost::counting_iterator<size_t>{input.size()},
                       output.begin(),
                       rw_tr);
    };

    struct frenel_inf {
        double out_w;
        double out_h;
        double in_h;
        double in_w;
        double lambda;
        double d;
    };

    
    std::valarray<std::complex<double>> fft(std::valarray<std::complex<double>>& input) {
        
        size_t n = input.size();
        
        

        
    };

    vec_complex fft2d(vec_complex& input, size_t w, size_t h) {
        vec_complex output(w*h);

        return output;
    }

    vec_complex frenel(vec_complex& input, frenel_inf inf, size_t w, size_t h) {
    
        vec_complex output(w*h);

        double out_w_ratio = inf.out_w / w ;
        double out_h_ratio = inf.out_h / h;
        double in_w_ratio = inf.in_w / w;
        double in_h_ratio = inf.in_h / h;
        double k = 2 * M_PI / inf.lambda;
        double lambda = inf.lambda;
        double d = inf.d; 
        
        auto sum_up ( [=](const size_t out_index) {
            double x_out = (out_index % w + 1 - w / 2.) * out_w_ratio;
            double y_out = (out_index / w + 1 - h / 2.) * out_h_ratio;
            return [=](const size_t in_index, const complex c) -> complex {
                double x_in = (in_index % w + 1 - w / 2.) * in_w_ratio;
                double y_in = (in_index / w + 1 - h / 2.) * in_h_ratio;
                return c * std::exp(1i * M_PI * (x_in*x_in + y_in*y_in) / lambda / d)
                         * std::exp(-1i * 2. * M_PI * (x_out / w + y_out / h));    
            };

        } );    

        auto frenel_tr ( [=](const size_t out_index) -> complex {
            
            complex result = std::transform_reduce(boost::counting_iterator<size_t>{0},
                                         boost::counting_iterator<size_t>{w*h},    
                                         input.begin(),
                                         complex{},
                                         std::plus<complex>(),
                                         sum_up(out_index));

            double coef = lambda * d;
            double x_out = (out_index % w + 1 - w / 2.) * out_w_ratio;
            double y_out = (out_index / w + 1 - h / 2.) * out_h_ratio;

            return result  * std::exp(1i * 2. * M_PI * d / lambda)
                           * std::exp(1i * M_PI * (x_out * x_out + y_out * y_out) / coef)
                           / 1i / coef;
        });

        std::transform(boost::counting_iterator<size_t>{0},
                       boost::counting_iterator<size_t>{w*h},
                       output.begin(),
                       frenel_tr);
        
        return output;
   }
}

#endif // OPTICS_H
