#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include<string>

class Polynom {
public:
	Polynom(int degree) : coeff_(degree + 1, 0), degree_(degree) {
	}

	Polynom(const std::vector<int>& coeff) : coeff_(coeff), degree_(coeff.size() - 1) {
	}

	Polynom(const Polynom& rhs) : degree_(rhs.degree_), coeff_(rhs.coeff_) {
	}

	Polynom(Polynom&& rhs) : degree_(rhs.degree_), coeff_(std::move(rhs.coeff_)) {
		rhs.degree_ = 0;
	}

	Polynom& operator=(const Polynom& rhs) {
		this->degree_ = rhs.degree_;
		this->coeff_ = rhs.coeff_;

		return *this;
	}

	static Polynom get_add_neutral(int degree) {
		return Polynom(degree);
	}

	static Polynom get_mul_neutral(int degree) {
		Polynom result = Polynom(degree);
		result.coeff_[degree] = 1;
		return result;
	}

	Polynom operator+(const Polynom& rhs) const {
		Polynom result(std::max(rhs.degree_, this->degree_));
		auto coeff_first = rhs.coeff_;
		auto coeff_second = this->coeff_;

		std::reverse(coeff_first.begin(), coeff_first.end());
		std::reverse(coeff_second.begin(), coeff_second.end());

		for (int i = 0; i < std::max(rhs.degree_, this->degree_) + 1; ++i) {
			int operand_first = 0, operand_second = 0;
			if (i < rhs.degree_ + 1) {
				operand_first = coeff_first[i];
			}
			if (i < this->degree_ + 1) {
				operand_second = coeff_second[i];
			}
			result.coeff_[i] = operand_first ^ operand_second;
		}

		std::reverse(result.coeff_.begin(), result.coeff_.end());

		return result;
	}


	Polynom operator*(const Polynom& rhs) const {
		Polynom result(rhs.degree_ + this->degree_);
		for (int idx_first = 0; idx_first <= rhs.degree_; ++idx_first) {
			for (int idx_second = 0; idx_second <= this->degree_; ++idx_second) {
				result.coeff_[idx_first + idx_second] ^= rhs.coeff_[idx_first] * this->coeff_[idx_second];
			}
		}

		return result.restrict_polynom();
	}


	static Polynom mul(const Polynom& lhs, const Polynom& rhs, const Polynom& gen) {
		Polynom result(rhs.degree_ + lhs.degree_);
		for (int idx_first = 0; idx_first <= rhs.degree_; ++idx_first) {
			for (int idx_second = 0; idx_second <= lhs.degree_; ++idx_second) {
				result.coeff_[idx_first + idx_second] ^= rhs.coeff_[idx_first] * lhs.coeff_[idx_second];
			}
		}

		return Polynom::reduction(result, gen).restrict_polynom();
	}

	static Polynom square(const Polynom& poly, const Polynom& gen) {
		return mul(poly, poly, gen);
	}

	static Polynom pow(const Polynom& poly, int power, const Polynom& gen) {
		if (power == 0) {
			return Polynom::polynom_of_degree(0);
		}
		else if (power == 1) {
			return poly;
		}
		else if (power % 2 == 0) {
			auto result = pow(poly, power / 2, gen);
			return Polynom::mul(result, result, gen);
		}
		else {
			return Polynom::mul(Polynom::pow(poly, power - 1, gen), poly, gen);
		}
	}

	static Polynom pow2(const Polynom& poly, int deg_two, const Polynom& gen) {
		if (deg_two == 0) {
			return poly;
		}
		else {
			auto result = pow2(poly, deg_two - 1, gen);
			return Polynom::mul(result, result, gen);
		}
	}

	static Polynom reverse_mul(const Polynom& poly, int m, const Polynom& gen) {
		Polynom result = Polynom::get_mul_neutral(poly.get_degree());
		Polynom current = Polynom::mul(poly, poly, gen);
		for (int i = m - 1; i >= 1; --i) {
			result = Polynom::mul(result, current, gen);
			current = Polynom::mul(current, current, gen);
			//result = Polynom::mul(result, pow2(poly, i, gen), gen);
		}

		return result;
	}

	static int get_ld(const Polynom& poly) {
		for (int i = 0; i < poly.coeff_.size(); ++i) {
			if (poly.coeff_[i]) {
				return poly.degree_ - i;
			}
		}
		return 0;
	}

	static Polynom polynom_of_degree(int deg) {
		Polynom result(deg);
		result.coeff_[0] = 1;
		return result;
	}

	static Polynom initialize_from_list_of_degree(const std::vector<int>& degrees) {
		int max_deg = *std::max_element(degrees.begin(), degrees.end());
		Polynom result(max_deg);
		for (const auto& deg : degrees) {
			result.coeff_[max_deg - deg] = 1;
		}

		return result;
	}

	static Polynom reduction(const Polynom& poly, const Polynom& mod) {
		Polynom result(poly);
		while (Polynom::get_ld(result) >= Polynom::get_ld(mod)) {
			int mul_degree = get_ld(result) - get_ld(mod);
			Polynom multiplier = Polynom::polynom_of_degree(mul_degree);
			result = result + (mod * multiplier);
		}
		return result;
	}

	void print() const {
		int last_one_idx = -1;
		for (int i = 0; i < coeff_.size(); ++i) {
			if (coeff_[i]) {
				last_one_idx = i;
			}
		}
		if (last_one_idx == -1) {
			std::cout << 0 << std::endl;
			return;
		}
		for (int i = 0; i < coeff_.size(); ++i) {
			if (coeff_[i] && i == static_cast<int>(coeff_.size()) - 1) {
				std::cout << "1";
			}
			else if (coeff_[i] && i == static_cast<int>(coeff_.size()) - 2) {
				std::cout << "x";
			}
			else if (coeff_[i]) {
				std::cout << "x^" << (static_cast<int>(coeff_.size()) - 1 - i);
			}
			if (coeff_[i] && i < last_one_idx) {
				std::cout << " + ";
			}
		}
		std::cout << std::endl;
	}

	int get_degree() const {
		return degree_;
	}

	std::string to_string(int m) const {
		std::string result;
		for (int i = 0; i < m - degree_ - 1; ++i) {
			result += "0";
		}
		for (int i = 0; i < coeff_.size(); ++i) {
			result += std::to_string(coeff_[i]);
		}

		return result;
	}

private:
	Polynom restrict_polynom() {
		bool flag = true;
		while (flag) {
			if (coeff_[0] == 0) {
				coeff_.erase(coeff_.begin());
				--degree_;
			}
			else {
				flag = false;
			}
		}

		return *this;
	}

	std::vector<int> coeff_;
	int degree_;
};



int main() {
	/*
		{
			const int m = 3;
			Polynom mod = Polynom({1, 0, 1, 1});
			mod.print();
			Polynom a = Polynom({1, 0, 0});

			a.print();
			auto rev = Polynom::reverse_mul(a, 3, mod);
			rev.print();
			Polynom::mul(a, rev, mod).print();
			return 0;
		}
	*/

	const int m = 509;
	Polynom generator = Polynom::initialize_from_list_of_degree({ 509, 23, 3, 2, 0 });
	Polynom a = Polynom::initialize_from_list_of_degree({ 6, 4, 2, 1, 0 });
	Polynom b = Polynom::initialize_from_list_of_degree({ 8, 7, 3, 1, 0 });
	Polynom c = Polynom::initialize_from_list_of_degree({ 7, 5, 2, 0 });

	std::cout << "generator = ";
	generator.print();
	std::cout << std::endl;

	std::cout << "a = ";
	a.print();
	std::cout << std::endl;

	std::cout << "b = ";
	b.print();
	std::cout << std::endl;

	std::cout << "c = ";
	c.print();
	std::cout << std::endl;

	std::cout << "a + b = ";
	(a + b).print();
	std::cout << std::endl;

	std::cout << "a * b = ";
	Polynom::mul(a, b, generator).print();
	std::cout << std::endl;

	std::cout << "a^2 = ";
	Polynom::square(a, generator).print();
	std::cout << std::endl;

	std::cout << "b^2 = ";
	Polynom::square(b, generator).print();
	std::cout << std::endl;

	std::cout << "a^25 = ";
	Polynom::pow(a, 25, generator).print();
	std::cout << std::endl;

	std::cout << "b^25 = ";
	Polynom::pow(b, 25, generator).print();
	std::cout << std::endl;

	std::cout << "a^(-1) = ";
	auto reverse_a = Polynom::reverse_mul(a, m, generator);
	reverse_a.print();
	std::cout << std::endl;

	std::cout << "a^(-1) * a = ";
	Polynom::mul(reverse_a, a, generator).print();
	std::cout << std::endl;

	std::cout << "b^(-1) = ";
	Polynom::reverse_mul(b, m, generator).print();
	std::cout << std::endl;

	std::cout << "b^(-1) * b = ";
	Polynom::mul(b, Polynom::reverse_mul(b, m, generator), generator).print();
	std::cout << std::endl;

	std::cout << "a = ";
	std::cout << a.to_string(m) << std::endl;
	std::cout << std::endl;

	std::cout << "b = ";
	std::cout << b.to_string(m) << std::endl;
	std::cout << std::endl;

	std::cout << "(a + b) * c = ";
	Polynom::mul(a + b, c, generator).print();
	std::cout << std::endl;

	std::cout << "a * c + b * c = ";
	(Polynom::mul(a, c, generator) + Polynom::mul(b, c, generator)).print();
	return 0;
}