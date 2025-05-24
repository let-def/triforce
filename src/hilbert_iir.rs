const ORDER: usize = 12;

#[derive(Clone, Copy)]
pub struct Filter {
    coeffs_r: [f32; ORDER],
    coeffs_i: [f32; ORDER],
    poles_r: [f32; ORDER],
    poles_i: [f32; ORDER],
    direct: f32,
}

#[derive(Clone, Copy)]
pub struct State {
    real: [f32; ORDER],
    imag: [f32; ORDER],
}

impl Filter {
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

    fn process<const N: usize>(state: &mut [State; N], input: &[&[f32]; N], output: &[&mut [Complex<f32>]]) {
    }
}
