/// Calculates value for B-spline based on knots, coefficients, order and position
pub fn deboor_alg(x: f64, knots: &[f64], coeffs: &[f64], order: &usize) -> f64 {
    // Assume that knots are sorted into increasing order
    let k: usize = (1..knots.len())
        .find(|&i| knots[i] >= x)
        .expect("Evaluation point does not lie within knot intervals")
        + 1;

    let p: usize = order - 1;
    let mut i: usize;
    let mut alpha_j: f64;
    let mut d_coeff: Vec<f64> = (0..=p).map(|i| coeffs[i + k - p]).collect();

    for r in 1..=p {
        for j in (r..=p).rev() {
            i = j + k - p;
            alpha_j = (x - knots[i]) / (knots[j + 1 + k - r] - knots[i]);
            d_coeff[j] = (1.0 - alpha_j) * d_coeff[j - 1] + alpha_j * d_coeff[j];
        }
    }

    d_coeff[p]
}
