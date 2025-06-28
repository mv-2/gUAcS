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

/// 2x2 matrix struct to simplify beam tracing formulae
#[derive(Clone, Copy)]
pub struct Mat2<T> {
    pub a: T,
    pub b: T,
    pub c: T,
    pub d: T,
}

impl std::ops::Add for Mat2<f64> {
    type Output = Mat2<f64>;

    fn add(self, _rhs: Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self.a + _rhs.a,
            b: self.b + _rhs.b,
            c: self.c + _rhs.c,
            d: self.d + _rhs.d,
        }
    }
}

impl std::ops::Sub for Mat2<f64> {
    type Output = Mat2<f64>;

    fn sub(self, _rhs: Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self.a - _rhs.a,
            b: self.b - _rhs.b,
            c: self.c - _rhs.c,
            d: self.d - _rhs.d,
        }
    }
}

impl std::ops::Mul for Mat2<f64> {
    type Output = Mat2<f64>;

    fn mul(self, _rhs: Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self.a * _rhs.a + self.b * _rhs.c,
            b: self.a * _rhs.b + self.b * _rhs.d,
            c: self.c * _rhs.a + self.a * _rhs.c,
            d: self.c * _rhs.b + self.d * _rhs.d,
        }
    }
}

impl std::ops::Mul<Mat2<f64>> for f64 {
    type Output = Mat2<f64>;

    fn mul(self, _rhs: Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self * _rhs.a,
            b: self * _rhs.b,
            c: self * _rhs.c,
            d: self * _rhs.d,
        }
    }
}

impl std::ops::Div<f64> for Mat2<f64> {
    type Output = Mat2<f64>;

    fn div(self, _rhs: f64) -> Mat2<f64> {
        Mat2 {
            a: self.a / _rhs,
            b: self.b / _rhs,
            c: self.c / _rhs,
            d: self.d / _rhs,
        }
    }
}

impl Mat2<f64> {
    #[inline]
    pub const fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Mat2 { a, b, c, d }
    }

    pub const I: Self = Self::new(1_f64, 0_f64, 0_f64, 1_f64);

    /// returns inverted [`Mat2`] struct
    pub fn inv(&self) -> Mat2<f64> {
        let det = (self.a * self.d - self.b * self.c).powi(-1);
        Mat2 {
            a: self.d / det,
            b: -self.b / det,
            c: -self.c / det,
            d: self.a / det,
        }
    }

    /// Matrix multiplication between self and second [`Mat2`]
    pub fn mul22(&self, mulmat: &Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self.a * mulmat.a + self.b * mulmat.c,
            b: self.a * mulmat.b + self.b * mulmat.d,
            c: self.c * mulmat.a + self.d * mulmat.c,
            d: self.c * mulmat.b + self.d * mulmat.d,
        }
    }

    pub fn add22(&self, addmat: &Mat2<f64>) -> Mat2<f64> {
        Mat2 {
            a: self.a + addmat.a,
            b: self.b + addmat.b,
            c: self.c + addmat.c,
            d: self.d + addmat.d,
        }
    }

    /// displays matrix
    pub fn disp(&self) {
        println!("[{:}, {:};\n {:}, {:}]", self.a, self.b, self.c, self.d);
    }
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

    #[test]
    fn matmul_test() {
        let ident: Mat2<f64> = Mat2::I;
        let prod_ident: Mat2<f64> = ident * ident;
        prod_ident.disp();
        assert!(
            (prod_ident.a == 1_f64)
                && (prod_ident.b == 0_f64)
                && (prod_ident.c == 0_f64)
                && (prod_ident.d == 1_f64)
        );
    }

    #[test]
    fn matadd_test() {
        let ident: Mat2<f64> = Mat2::I;
        let add_ident: Mat2<f64> = ident + ident;
        add_ident.disp();
        assert!(
            (add_ident.a == 2_f64)
                && (add_ident.b == 0_f64)
                && (add_ident.c == 0_f64)
                && (add_ident.d == 2_f64)
        );
    }

    #[test]
    fn matsub_test() {
        let ident: Mat2<f64> = Mat2::I;
        let sub_ident: Mat2<f64> = ident - ident;
        sub_ident.disp();
        assert!(
            (sub_ident.a == 0_f64)
                && (sub_ident.b == 0_f64)
                && (sub_ident.c == 0_f64)
                && (sub_ident.d == 0_f64)
        );
    }

    #[test]
    fn mat_scalar_prod_test() {
        let ident: Mat2<f64> = Mat2::I;
        let scal_prod: Mat2<f64> = 5_f64 * ident;
        scal_prod.disp();
        assert!(
            (scal_prod.a == 5_f64)
                && (scal_prod.b == 0_f64)
                && (scal_prod.c == 0_f64)
                && (scal_prod.d == 5_f64)
        );
    }

    #[test]
    fn mat_scalar_div_test() {
        let ident: Mat2<f64> = Mat2::I;
        let scal_div: Mat2<f64> = ident / 2_f64;
        scal_div.disp();
        assert!(
            (scal_div.a == 0.5)
                && (scal_div.b == 0_f64)
                && (scal_div.c == 0_f64)
                && (scal_div.d == 0.5)
        );
    }
}
