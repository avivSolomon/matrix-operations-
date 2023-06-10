
class Polynomial:
    def __init__(self, coeffs):
        self.coeffs = coeffs

    def evaluate(self, x):
        """
        Evaluates the polynomial at a given value of x.
        """
        result = 0
        for i in range(len(self.coeffs)):
            result += self.coeffs[i] * x**i
        return result

    def derivative(self):
        """
        Computes the derivative of the polynomial.
        """
        derivative_coeffs = [i * self.coeffs[i] for i in range(1, len(self.coeffs))]
        return Polynomial(derivative_coeffs)

    def newton_raphson(self, initial_guess, epsilon=1e-6, max_iterations=100):
        """
        Performs the Newton-Raphson method to find the roots of the polynomial.

        Args:
            initial_guess (float): The initial guess for the root.
            epsilon (float): The desired accuracy of the root. Defaults to 1e-6.
            max_iterations (int): The maximum number of iterations. Defaults to 100.

        Returns:
            float: The estimated root of the polynomial.
        """
        x = initial_guess
        iterations = 0

        while abs(self.evaluate(x)) > epsilon and iterations < max_iterations:
            x = x - self.evaluate(x) / self.derivative().evaluate(x)
            iterations += 1
        return x