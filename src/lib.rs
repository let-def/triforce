// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * An attempt at an MVDR beamformer for equilateral triangle mic arrays
 * as found on newer Apple Silicon Macs.
 *
 * Currently mono, but could probably be extended to be stereo somewhat
 * easily. I think.
 *
 * Copyright (C) 2024 James Calligeros <jcalligeros99@gmail.com>
 */

use std::f32::consts::PI;
use itertools::izip;
use lv2::prelude::*;

use nalgebra::{linalg::SVD, Complex, Matrix3, Vector3};

const C: f32 = 343.00; /* m*s^-1 */

const ORDER: usize = 12;

#[derive(Clone, Copy)]
pub struct Coeffs {
    coeffs_r: [f32; ORDER],
    coeffs_i: [f32; ORDER],
    poles_r: [f32; ORDER],
    poles_i: [f32; ORDER],
    direct: f32,
}

#[derive(Clone, Copy)]
pub struct Filter {
    real: [f32; ORDER],
    imag: [f32; ORDER],
}

impl Coeffs {
    fn init(sample_rate: f32, passband_gain: f32) -> Coeffs {
        let freq_factor = f32::min(0.46, 20000.0 / sample_rate);
        let mut result = Coeffs {
            coeffs_r: [0.0; ORDER],
            coeffs_i: [0.0; ORDER],
            poles_r: [0.0; ORDER],
            poles_i: [0.0; ORDER],
            direct: 0.0,
        };
        result.direct = CONSTANTS.direct * 2.0 * passband_gain * freq_factor;
        for i in 0..ORDER {
            let coeff = CONSTANTS.coeffs[i] * freq_factor * passband_gain;
            result.coeffs_r[i] = coeff.re;
            result.coeffs_i[i] = coeff.im;
            let pole = (CONSTANTS.poles[i] * freq_factor).exp();
            result.poles_r[i] = pole.re;
            result.poles_i[i] = pole.im;
        }
        result
    }
}

impl Filter {
    const INITIAL : Filter = Filter {
        real: [-1.0; ORDER],
        imag: [0.0; ORDER],
    };

    // fn sample(&mut self, coeffs: &Coeffs, x: f32) -> (f32, f32) {
    //     let mut result_r = x * coeffs.direct;
    //     let mut result_i = 0.0;
    //     for i in 0..ORDER {
    //         let real = self.real[i] * coeffs.poles_r[i] - self.imag[i] * coeffs.poles_i[i] + x * coeffs.coeffs_r[i];
    //         let imag = self.real[i] * coeffs.poles_i[i] + self.imag[i] * coeffs.poles_r[i] + x * coeffs.coeffs_i[i];
    //         result_r += real;
    //         self.real[i] = real;
    //         result_i += imag;
    //         self.imag[i] = imag;
    //     }
    //     (result_r, result_i)
    // }

    // fn process(&mut self, coeffs: &Coeffs, x: &[f32], r: &mut [f32], i: &mut [f32]) {
    //     for (x, r, i) in izip!(x, r, i) {
    //         let (rr, ri) = self.sample(coeffs, *x);
    //         *r = rr;
    //         *i = ri;
    //     }
    // }

    // fn process2(
    //     coeffs: &Coeffs,
    //     state: &mut [Self; 2],
    //     x0: &[f32],
    //     x1: &[f32],
    //     r0: &mut [f32],
    //     r1: &mut [f32],
    //     i0: &mut [f32],
    //     i1: &mut [f32],
    // ) {
    //     for (x0, x1, r0, r1, i0, i1) in izip!(x0, x1, r0, r1, i0, i1) {
    //         let mut r0v = x0 * coeffs.direct;
    //         let mut i0v = 0.0;
    //         let mut r1v = x1 * coeffs.direct;
    //         let mut i1v = 0.0;
    //         for j in 0..ORDER {
    //             let real0 = state[0].real[j] * coeffs.poles_r[j]
    //                 - state[0].imag[j] * coeffs.poles_i[j]
    //                 + x0 * coeffs.coeffs_r[j];
    //             let imag0 = state[0].real[j] * coeffs.poles_i[j]
    //                 + state[0].imag[j] * coeffs.poles_r[j]
    //                 + x0 * coeffs.coeffs_i[j];
    //             let real1 = state[1].real[j] * coeffs.poles_r[j]
    //                 - state[1].imag[j] * coeffs.poles_i[j]
    //                 + x1 * coeffs.coeffs_r[j];
    //             let imag1 = state[1].real[j] * coeffs.poles_i[j]
    //                 + state[1].imag[j] * coeffs.poles_r[j]
    //                 + x1 * coeffs.coeffs_i[j];
    //             r0v += real0;
    //             state[0].real[j] = real0;
    //             i0v += imag0;
    //             state[0].imag[j] = imag0;
    //             r1v += real1;
    //             state[1].real[j] = real1;
    //             i1v += imag1;
    //             state[1].imag[j] = imag1;
    //         }
    //         *r0 = r0v;
    //         *i0 = i0v;
    //         *r1 = r1v;
    //         *i1 = i1v;
    //     }
    // }

    fn process3(
        coeffs: &Coeffs,
        state: &mut [Self; 3],
        x0: &[f32],
        x1: &[f32],
        x2: &[f32],
        r0: &mut [f32],
        r1: &mut [f32],
        r2: &mut [f32],
        i0: &mut [f32],
        i1: &mut [f32],
        i2: &mut [f32],
    ) {
        for (x0, x1, x2, r0, r1, r2, i0, i1, i2) in izip!(x0, x1, x2, r0, r1, r2, i0, i1, i2) {
            let mut r0_val = x0 * coeffs.direct;
            let mut i0_val = 0.0;
            let mut r1_val = x1 * coeffs.direct;
            let mut i1_val = 0.0;
            let mut r2_val = x2 * coeffs.direct;
            let mut i2_val = 0.0;
            for j in 0..ORDER {
                let real0 = state[0].real[j] * coeffs.poles_r[j]
                    - state[0].imag[j] * coeffs.poles_i[j]
                    + x0 * coeffs.coeffs_r[j];
                let imag0 = state[0].real[j] * coeffs.poles_i[j]
                    + state[0].imag[j] * coeffs.poles_r[j]
                    + x0 * coeffs.coeffs_i[j];
                let real1 = state[1].real[j] * coeffs.poles_r[j]
                    - state[1].imag[j] * coeffs.poles_i[j]
                    + x1 * coeffs.coeffs_r[j];
                let imag1 = state[1].real[j] * coeffs.poles_i[j]
                    + state[1].imag[j] * coeffs.poles_r[j]
                    + x1 * coeffs.coeffs_i[j];
                let real2 = state[2].real[j] * coeffs.poles_r[j]
                    - state[2].imag[j] * coeffs.poles_i[j]
                    + x2 * coeffs.coeffs_r[j];
                let imag2 = state[2].real[j] * coeffs.poles_i[j]
                    + state[2].imag[j] * coeffs.poles_r[j]
                    + x2 * coeffs.coeffs_i[j];
                r0_val += real0;
                state[0].real[j] = real0;
                i0_val += imag0;
                state[0].imag[j] = imag0;
                r1_val += real1;
                state[1].real[j] = real1;
                i1_val += imag1;
                state[1].imag[j] = imag1;
                r2_val += real2;
                state[2].real[j] = real2;
                i2_val += imag2;
                state[2].imag[j] = imag2;
            }
            *r0 = r0_val;
            *i0 = i0_val;
            *r1 = r1_val;
            *i1 = i1_val;
            *r2 = r2_val;
            *i2 = i2_val;
        }
    }

}

struct Constants {
    coeffs: [Complex<f32>; ORDER],
    poles: [Complex<f32>; ORDER],
    direct: f32,
}

const CONSTANTS: Constants = Constants {
    coeffs: [
        Complex::<f32>::new(-0.000224352093802, 0.00543499018201),
        Complex::<f32>::new(0.010750055781500, -0.01738906856810),
        Complex::<f32>::new(-0.045679587391700, 0.02291669314290),
        Complex::<f32>::new(0.112825005820000, 0.00278413661237),
        Complex::<f32>::new(-0.208067578452000, -0.10462895867500),
        Complex::<f32>::new(0.287178375010000, 0.33619239719000),
        Complex::<f32>::new(-0.254675294431000, -0.68303389965500),
        Complex::<f32>::new(0.048108183502600, 0.95406158937400),
        Complex::<f32>::new(0.227861357867000, -0.89127357456900),
        Complex::<f32>::new(-0.365411839137000, 0.52508831727100),
        Complex::<f32>::new(0.280729061131000, -0.15513120660600),
        Complex::<f32>::new(-0.093506178772800, 0.00512245855404),
    ],
    poles: [
        Complex::<f32>::new(-0.00495335976478, 0.0092579876872),
        Complex::<f32>::new(-0.01785949130200, 0.0273493725543),
        Complex::<f32>::new(-0.04137143731550, 0.0744756910287),
        Complex::<f32>::new(-0.08821484088850, 0.1783496774570),
        Complex::<f32>::new(-0.17922965812000, 0.3960134022300),
        Complex::<f32>::new(-0.33826180075300, 0.8292295333540),
        Complex::<f32>::new(-0.55768869973200, 1.6129853832800),
        Complex::<f32>::new(-0.73515773614800, 2.7998739868200),
        Complex::<f32>::new(-0.71905738117200, 4.1639616612800),
        Complex::<f32>::new(-0.51787102520900, 5.2972482680400),
        Complex::<f32>::new(-0.28019746947100, 5.9959860238800),
        Complex::<f32>::new(-0.08527513545310, 6.3048492377000),
    ],
    direct: 0.000262057212648,
};

/// The distance of a given element in the array from the zeroth
/// element
#[derive(Copy, Clone)]
struct ElemDistance {
    x: f32,
    y: f32,
}

/// The steering vector is a representation of the phase delays at each microphone.
/// It is calculated by taking the dot product of the array geometry matrix and the
/// unit vector of the direction of arrival.
fn steering_vec(theta: f32, phi: f32, f: f32, elems: [ElemDistance; 3]) -> Vector3<Complex<f32>> {
    // Mic positions are relative to Left/Top to preserve x/y axis semantics
    let mic_positions: Vec<Vector3<f32>> =
        elems.iter().map(|e| Vector3::new(e.x, e.y, 0f32)).collect();

    // Calculate angular repetency (2pi/lambda)
    let repetency = (2f32 * PI) / (f / C);

    // Compute the unit vector of the DOA
    let u_dir = Vector3::new(phi.sin() * theta.cos(), phi.sin() * theta.sin(), phi.cos());

    // Calculate the steering vector by taking the array geometry, speed of sound,
    // and DOA unit vector
    let mut steering_vector = Vector3::from_element(Complex::new(0f32, 0f32));

    for (i, mic_pos) in mic_positions.iter().enumerate() {
        let delay = mic_pos.dot(&u_dir) / C;
        let phase = -repetency * delay;
        steering_vector[i] = Complex::new(phase.cos(), phase.sin());
    }

    steering_vector
}

/// There's nothing special about this, it's just a covariance matrix. It is always
/// square.
fn covariance(real: &[Vec<f32>; 3], imag: &[Vec<f32>; 3]) -> Matrix3<Complex<f32>> {
    let n_samples = real[0].len();

    let mut covar = Matrix3::zeros();

    for t in 0..n_samples {
        let discrete: Vector3<Complex<f32>> = Vector3::new(
            Complex::new(real[0][t], imag[0][t]),
            Complex::new(real[1][t], imag[1][t]),
            Complex::new(real[2][t], imag[2][t])
        );
        covar += &discrete * discrete.adjoint();
    }

    // Our samples are shit, so we can't get a very nice covariance matrix.
    // Regularise the shit covariance matrix by introducing a constant value
    // across the identity
    let reg = Matrix3::identity().map(|x: f32| Complex::new(x * 1e-4f32, 0f32));

    covar /= Complex::new(n_samples as f32, 0f32);
    covar + reg
}

/// To calculate the weighting vector for the beamformer, we need to invert
/// the covariance matrix, multiply it by the steering vector, then divide that
/// by itself multiplied by the Hermitian transpose of the steering vector
/// w = (cov^-1 * sv) / (sv.adjoint() * (cov^-1 * sv)). Note that the denominator
/// is the same as the conjugate-linear dot product of the steering vector and
/// the numerator.
fn mvdr_weights(cov: &Matrix3<Complex<f32>>, sv: &Vector3<Complex<f32>>) -> Vector3<Complex<f32>> {
    // Since we have a numerically unstable covariance matrix, we can't take the
    // true inverse of it. Let's instead decompose it and take the pseudoinverse.
    let svd = SVD::new(cov.to_owned(), true, true);
    let r_inv = svd.pseudo_inverse(1e-4f32).unwrap();

    let num = r_inv * sv;
    let den = sv.dotc(&num); // Conjugate-linear dot product

    num / den
}

/*
 * Input and output ports used by the plugin
 *
 *
 * Ports:
 *      in_1: channel 1 input (left/top)
 *      in_2: channel 2 input (right/bottom)
 *      in_3: channel 3 input (vertex)
 *      out: output
 *      h_angle: horizontal steering angle in degrees (relative to input 1)
 *      v_angle: vertical steering angle in degrees
 *      opt_freq: frequency to optimise for
 *      t_win: covariance matrix time window
 */
#[derive(PortCollection)]
pub struct Ports {
    in_1: InputPort<Audio>,
    in_2: InputPort<Audio>,
    in_3: InputPort<Audio>,
    out: OutputPort<Audio>,
    h_angle: InputPort<Control>,
    v_angle: InputPort<Control>,
    opt_freq: InputPort<Control>,
    t_win: InputPort<Control>,
    mic2_x: InputPort<Control>,
    mic2_y: InputPort<Control>,
    mic3_x: InputPort<Control>,
    mic3_y: InputPort<Control>,
}

/*
 * Plugin state
 */
#[uri("https://chadmed.au/triforce")]
pub struct Triforce {
    hangle_curr: f32,
    vangle_curr: f32,
    freq_curr: f32,
    sample_rate: f32,
    samples_since_last_update: usize,
    steering_vector: Vector3<Complex<f32>>,
    window_real: [Vec<f32>; 3],
    window_imag: [Vec<f32>; 3],
    covar: Matrix3<Complex<f32>>,
    array_geom: [ElemDistance; 3],
    weights: Vector3<Complex<f32>>,
    analytic_coeffs: Coeffs,
    analytic_filters: [Filter; 3],
    input_real0: Vec<f32>,
    input_real1: Vec<f32>,
    input_real2: Vec<f32>,
    input_imag0: Vec<f32>,
    input_imag1: Vec<f32>,
    input_imag2: Vec<f32>,
}

trait Beamformer: Plugin {
    fn update_params(&mut self, ports: &mut Ports);
}

impl Triforce {
    pub fn with_sample_rate(sample_rate: f32) -> Self {
        Self {
            hangle_curr: 0f32,
            vangle_curr: 0f32,
            freq_curr: 1000f32,
            samples_since_last_update: usize::max_value(),
            sample_rate,
            window_real: [Vec::new(),Vec::new(),Vec::new()],
            window_imag: [Vec::new(),Vec::new(),Vec::new()],
            array_geom: [ElemDistance { x: 0f32, y: 0f32 }; 3],
            steering_vector: steering_vec(
                90f32.to_radians(),
                45f32.to_radians(),
                1000f32,
                [ElemDistance { x: 0f32, y: 0f32 }; 3],
            ),
            covar: Matrix3::zeros(),
            weights: Vector3::zeros(),
            analytic_coeffs: Coeffs::init(sample_rate, 2.0),
            analytic_filters: [Filter::INITIAL,Filter::INITIAL,Filter::INITIAL],
            input_real0: Vec::new(),
            input_real1: Vec::new(),
            input_real2: Vec::new(),
            input_imag0: Vec::new(),
            input_imag1: Vec::new(),
            input_imag2: Vec::new(),
        }
    }

    pub fn process_slice(
        &mut self,
        mic1: &[f32],
        mic2: &[f32],
        mic3: &[f32],
        output: &mut [f32],
        t_win: f32,
        buf_len: usize
    ) {
        // Steering vector is relative to Left/Top mic
        self.input_real0.resize(buf_len, 0.0);
        self.input_imag0.resize(buf_len, 0.0);
        self.input_real1.resize(buf_len, 0.0);
        self.input_imag1.resize(buf_len, 0.0);
        self.input_real2.resize(buf_len, 0.0);
        self.input_imag2.resize(buf_len, 0.0);
        Filter::process3(&self.analytic_coeffs, &mut self.analytic_filters,
                         mic1, mic2, mic3,
                         &mut self.input_real0,
                         &mut self.input_real1,
                         &mut self.input_real2,
                         &mut self.input_imag0,
                         &mut self.input_imag1,
                         &mut self.input_imag2);

        // Update the covariance matrix. We use an overlapping window to smooth over
        // the transitions.
        if self.samples_since_last_update as f32 >= (t_win / 1000f32) * self.sample_rate {
            self.samples_since_last_update = 0;
            // We want a 1/3 overlap
            let i = buf_len / 3;
            self.window_real[0].extend_from_slice(&self.input_real0[0..i]);
            self.window_real[1].extend_from_slice(&self.input_real1[0..i]);
            self.window_real[2].extend_from_slice(&self.input_real2[0..i]);
            self.window_imag[0].extend_from_slice(&self.input_imag0[0..i]);
            self.window_imag[1].extend_from_slice(&self.input_imag1[0..i]);
            self.window_imag[2].extend_from_slice(&self.input_imag2[0..i]);
            self.covar = covariance(&self.window_real, &self.window_imag);
            self.window_real[0] = self.input_real0[i..buf_len].to_vec();
            self.window_real[1] = self.input_real1[i..buf_len].to_vec();
            self.window_real[2] = self.input_real2[i..buf_len].to_vec();
            self.window_imag[0] = self.input_imag0[i..buf_len].to_vec();
            self.window_imag[1] = self.input_imag1[i..buf_len].to_vec();
            self.window_imag[2] = self.input_imag2[i..buf_len].to_vec();
            self.weights = mvdr_weights(&self.covar, &self.steering_vector);
        } else {
            self.samples_since_last_update += buf_len;
        }

        for t in 0..buf_len {
            let discrete: Vector3<Complex<f32>> = Vector3::new(
                Complex::new(self.input_real0[t], self.input_imag0[t]),
                Complex::new(self.input_real1[t], self.input_imag1[t]),
                Complex::new(self.input_real2[t], self.input_imag2[t])
            );

            let out =
                // Conjugate-linear dot product
                self.weights.dotc(&discrete)
                // Now we need to revert the Hilbert transform and output the signal
                .re;

            // Do all of our NFP and clamping here
            output[t] = if out.is_finite() && !out.is_nan() {
                out.clamp(-10f32, 10f32)
            } else {
                0f32
            };
        }
    }
}

impl Plugin for Triforce {
    type Ports = Ports;

    type InitFeatures = ();
    type AudioFeatures = ();

    fn new(info: &PluginInfo, _features: &mut ()) -> Option<Self> {
        Some(Self::with_sample_rate(info.sample_rate() as f32))
    }

    fn run(&mut self, ports: &mut Ports, _features: &mut (), _: u32) {
        Beamformer::update_params(self, ports);
        self.process_slice(
            &ports.in_1,
            &ports.in_2,
            &ports.in_3,
            &mut ports.out,
            *ports.t_win,
            samples as usize
        );
    }
}

impl Beamformer for Triforce {
    fn update_params(&mut self, ports: &mut Ports) {
        if self.hangle_curr != *ports.h_angle
            || self.freq_curr != *ports.opt_freq
            || self.vangle_curr != *ports.v_angle
            || self.array_geom[1].x != *ports.mic2_x
            || self.array_geom[1].y != *ports.mic2_y
            || self.array_geom[2].x != *ports.mic3_x
            || self.array_geom[2].y != *ports.mic3_y
        {
            self.hangle_curr = *ports.h_angle;
            self.vangle_curr = *ports.v_angle;
            self.freq_curr = *ports.opt_freq;
            self.array_geom = [
                ElemDistance { x: 0f32, y: 0f32 },
                ElemDistance {
                    x: *ports.mic2_x,
                    y: *ports.mic2_y,
                },
                ElemDistance {
                    x: *ports.mic3_x,
                    y: *ports.mic3_y,
                },
            ];
            self.steering_vector = steering_vec(
                self.hangle_curr.to_radians(),
                self.vangle_curr.to_radians(),
                self.freq_curr,
                self.array_geom,
            );

            // The steering vector has changed
            self.weights = mvdr_weights(&self.covar, &self.steering_vector);
        }
    }
}

lv2_descriptors!(Triforce);
