/// Calculates value for B-spline based on knots, coefficients, order and position
pub fn deboor_alg(x: f64, knots: &[f64], coeffs: &[f64], order: usize) -> f64 {
    // Assume that knots are sorted into increasing order
    let k: usize = (1..knots.len())
        .find(|&i| knots[i] >= x)
        .expect("Evaluation point does not lie within knot intervals")
        - 1;

    let mut i: usize;
    let mut alpha_j: f64;
    let mut d_coeff: Vec<f64> = (0..=order).map(|j| coeffs[j + k - order]).collect();

    for r in 1..=order {
        for j in (r..=order).rev() {
            i = j + k - order;
            alpha_j = (x - knots[i]) / (knots[j + 1 + k - r] - knots[i]);
            d_coeff[j] = (1.0 - alpha_j) * d_coeff[j - 1] + alpha_j * d_coeff[j];
        }
    }

    d_coeff[order]
}

/// RK4 step for variable step size
pub fn rk4_var_step_2d() -> [f64; 2] {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn deboor_test() {
        const POINTS: [f64; 4] = [1.5, 30.0, 57.0, 99.0];
        const VALS: [f64; 4] = [-6.3875, -104.0, -121.55, -3.95];
        let knots: Vec<f64> = vec![
            0.0, 0.0, 0.0, 0.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 100.0, 100.0,
            100.0,
        ];
        let coefs: Vec<f64> = vec![
            1.0,
            -32.33333333,
            -72.33333333,
            -105.66666667,
            -120.66666667,
            -125.66666667,
            -120.66666667,
            -105.66666667,
            -72.33333333,
            -32.33333333,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let order: usize = 3;
        let calc_vals: Vec<f64> = POINTS
            .iter()
            .map(|&x| deboor_alg(x, &knots, &coefs, order))
            .collect::<Vec<f64>>();
        // println!("{:#?}", calc_vals);
        for (i, &val) in calc_vals.iter().enumerate() {
            println!("{val}, {:}", VALS[i]);
            assert!((val - VALS[i]).abs() < 1e-8);
        }
    }
}
