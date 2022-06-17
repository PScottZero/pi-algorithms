using Plots
using Printf

function pi_chudnovsky(num_terms)
  sum_of_terms = sum([chudnovsky_term(big(q - 1)) for q in 1:num_terms])
  return (426880 * sqrt(big(10005))) / sum_of_terms
end

function chudnovsky_term(q)
  numerator = factorial(6 * q) * (545140134 * q + 13591409)
  denominator = factorial(3 * q) * (factorial(q) ^ 3) * ((-262537412640768000) ^ q)
  return numerator / denominator
end

function chudnovsky_iteration(pi_estimate, q)
  numerator = factorial(6 * q) * (545140134 * q + 13591409)
  denominator = factorial(3 * q) * (factorial(q) ^ 3) * ((-262537412640768000) ^ q)
  return pi_estimate + (426880 * sqrt(big(10005))) / (numerator / denominator)
end

function pi_gauss_legendre(num_iterations)
  a = big(1)
  b = big(1) / sqrt(big(2))
  t = big(1) / big(4)
  p = big(1)
  for _ in 1:num_iterations
    _, a, b, t, p = gauss_legendre_iteration(a, b, t, p)
  end
  return ((a + b) ^ 2) / (4 * t)
end

function gauss_legendre_iteration(a, b, t, p)
  a1 = (a + b) / 2
  b = sqrt(a * b)
  t = t - p * (a - a1) ^ 2
  p = 2 * p
  a = a1
  return ((a + b) ^ 2) / (4 * t), a, b, t, p
end

function digits_correct(pi_estimate)
  correct = 0
  pi_actual = replace(readline(open("pi-billion.txt", "r")), "." => "")
  pi_estimate = rpad(replace(string(pi_estimate), "." => ""), length(pi_actual), "0")
  for i in 1:length(pi_actual)
    if pi_estimate[i] != pi_actual[i]
      break
    end
    correct += 1
  end
  return correct
end

setprecision(4000000000)
pi_estimate = @time pi_gauss_legendre(30)
correct = digits_correct(pi_estimate)
println(correct)
