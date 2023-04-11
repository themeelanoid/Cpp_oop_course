#include <fstream>
#include <iostream>
#include <string>
#include <vector>

enum class Sign { PLUS, MINUS };

class BigInteger {
  friend std::istream& operator>>(std::istream& in, BigInteger& big_int);
  friend std::ostream& operator<<(std::ostream& out, const BigInteger& big_int);
  friend bool operator==(const BigInteger& left, const BigInteger& right);
  friend bool operator<(const BigInteger& left, const BigInteger& right);

 public:
  BigInteger();
  BigInteger(int number);
  BigInteger& operator=(const BigInteger& other);
  std::string toString() const;
  BigInteger operator-() const;
  BigInteger& operator+=(const BigInteger& other);
  BigInteger& operator++();
  BigInteger operator++(int);
  BigInteger& operator-=(const BigInteger& other);
  BigInteger& operator--();
  BigInteger operator--(int);
  BigInteger& operator*=(const BigInteger& other);
  void shift();
  BigInteger& operator/=(const BigInteger& other);
  BigInteger& operator%=(const BigInteger& other);
  explicit operator bool() const;
  bool even();
  void changeSign();
  BigInteger abs() const;

 private:
  std::vector<int> digits;
  static const int kDigitLimit = 1000000000;
  void removeLeadingZeros();
  void normalizeSign();
  void exchange(BigInteger& other);
  Sign sign;
};

BigInteger factorial(BigInteger other) {
  BigInteger a = 0;
  for (int i = 1; i < other; ++i) {
    a += i;
  }
  return a;
}

void BigInteger::removeLeadingZeros() {
  if (!digits.empty()) {
    size_t i = digits.size() - 1;
    while (digits[i] == 0 && i > 0) {
      digits.pop_back();
      --i;
    }
    normalizeSign();
  }
}

void BigInteger::normalizeSign() {
  if (digits.size() == 1 && digits[0] == 0) sign = Sign::PLUS;
}

void BigInteger::exchange(BigInteger& other) {
  std::swap(sign, other.sign);
  std::swap(digits, other.digits);
}

bool BigInteger::even() { return digits[0] % 2 == 0; }

BigInteger& BigInteger::operator=(const BigInteger& other) {
  BigInteger copy = other;
  exchange(copy);
  return *this;
}

std::string BigInteger::toString() const {
  std::string str;
  if (sign == Sign::MINUS) str.push_back('-');
  for (int i = static_cast<int>(digits.size() - 1); i >= 0; --i) {
    std::string digit = std::to_string(digits[i]);
    if ((i != static_cast<int>(digits.size() - 1)) && digit.length() < 9) {
      for (size_t j = 0; j < (9 - digit.length()); ++j) {
        str.push_back('0');
      }
    }
    str += digit;
  }
  return str;
}

BigInteger::BigInteger() : sign(Sign::PLUS) {}

BigInteger::BigInteger(int number) {
  if (number == 0) digits.push_back(number);
  sign = (number >= 0) ? Sign::PLUS : Sign::MINUS;
  number = (number < 0) ? -number : number;
  while (number != 0) {
    digits.push_back(number % kDigitLimit);
    number /= kDigitLimit;
  }
}

std::istream& operator>>(std::istream& in, BigInteger& big_int) {
  big_int.sign = Sign::PLUS;
  big_int.digits.clear();
  std::string str;
  in >> str;
  if (str.empty())
    big_int.sign = Sign::PLUS;
  else {
    if (str[0] == '-') {
      big_int.sign = Sign::MINUS;
      str = str.substr(1, str.length() - 1);
    }

    for (int i = str.length(); i > 0; i -= 9) {
      if (i < 9)
        big_int.digits.push_back(std::stoi(str.substr(0, i)));
      else
        big_int.digits.push_back(std::stoi(str.substr(i - 9, 9)));
    }
  }
  big_int.removeLeadingZeros();
  return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& big_int) {
  out << big_int.toString();
  return out;
}

bool operator==(const BigInteger& left, const BigInteger& right) {
  if (left.sign != right.sign || left.digits.size() != right.digits.size())
    return false;
  for (size_t i = 0; i < left.digits.size(); ++i) {
    if (left.digits[i] != right.digits[i]) {
      return false;
    }
  }
  return true;
}

bool operator<(const BigInteger& left, const BigInteger& right) {
  if (left == right) return false;
  if (left.sign == Sign::MINUS) {
    if (right.sign == Sign::PLUS) return true;
    return -right < -left;
  }
  if (right.sign == Sign::MINUS) return false;
  if (left.digits.size() > right.digits.size()) return false;
  if (left.digits.size() < right.digits.size()) return true;
  for (int i = static_cast<int>(left.digits.size()) - 1; i >= 0; --i) {
    if (left.digits[i] < right.digits[i]) return true;
    if (left.digits[i] > right.digits[i]) return false;
  }
  return false;
}

bool operator!=(const BigInteger& left, const BigInteger& right) {
  return !(left == right);
}

bool operator<=(const BigInteger& left, const BigInteger& right) {
  return !(right < left);
}

bool operator>(const BigInteger& left, const BigInteger& right) {
  return right < left;
}

bool operator>=(const BigInteger& left, const BigInteger& right) {
  return !(left < right);
}

void BigInteger::changeSign() {
  if (*this != 0) sign = sign == Sign::PLUS ? Sign::MINUS : Sign::PLUS;
}

BigInteger BigInteger::operator-() const {
  BigInteger copy = BigInteger(*this);
  if (copy != 0) copy.sign = copy.sign == Sign::PLUS ? Sign::MINUS : Sign::PLUS;
  return copy;
}

BigInteger BigInteger::abs() const {
  BigInteger copy = *this;
  if (copy.sign == Sign::MINUS) copy.changeSign();
  return copy;
}

BigInteger& BigInteger::operator+=(const BigInteger& other) {
  int carry = 0;
  if (sign == other.sign) {
    for (size_t i = 0;
         i < std::max(digits.size(), other.digits.size()) || carry != 0; ++i) {
      if (i == digits.size()) digits.push_back(0);
      digits[i] += carry + (i < other.digits.size() ? other.digits[i] : 0);
      carry = digits[i] >= kDigitLimit;
      if (carry != 0) digits[i] -= kDigitLimit;
    }
  } else {
    if (abs() >= other.abs()) {
      for (size_t i = 0; i < other.digits.size() || carry != 0; ++i) {
        digits[i] -= carry + (i < other.digits.size() ? other.digits[i] : 0);
        carry = digits[i] < 0;
        if (carry != 0) digits[i] += kDigitLimit;
      }
    } else {
      changeSign();
      *this = -other += *this;
      changeSign();
    }
    removeLeadingZeros();
  }
  normalizeSign();
  return *this;
}

BigInteger operator+(const BigInteger& left, const BigInteger& right) {
  BigInteger copy = left;
  copy += right;
  return copy;
}

BigInteger& BigInteger::operator++() {
  *this += 1;
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger copy = *this;
  *this += 1;
  return copy;
}

BigInteger& BigInteger::operator-=(const BigInteger& other) {
  changeSign();
  *this += other;
  changeSign();
  return *this;
}

BigInteger operator-(const BigInteger& left, const BigInteger& right) {
  BigInteger copy = left;
  copy -= right;
  return copy;
}

BigInteger& BigInteger::operator--() {
  *this -= 1;
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger copy = *this;
  *this -= 1;
  return copy;
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
  if (other == 0 || *this == 0) return *this = 0;
  if (other == 1) return *this;
  if (other == -1) {
    changeSign();
    return *this;
  }
  if (*this == 1) return *this = other;
  if (*this == -1) {
    *this = other;
    changeSign();
    return *this;
  }
  sign = sign == other.sign ? Sign::PLUS : Sign::MINUS;
  BigInteger copy = *this;
  digits.assign(digits.size() + other.digits.size(), 0);
  for (size_t i = 0; i < copy.digits.size(); ++i) {
    int carry = 0;
    for (size_t j = 0; j < other.digits.size() || carry != 0; ++j) {
      long long cur = digits[i + j] +
                      copy.digits[i] * 1LL *
                          (j < other.digits.size() ? other.digits[j] : 0) +
                      carry;
      digits[i + j] = static_cast<int>(cur % kDigitLimit);
      carry = static_cast<int>(cur / kDigitLimit);
    }
  }
  removeLeadingZeros();
  normalizeSign();
  return *this;
}

void BigInteger::shift() {
  if (digits.empty()) {
    digits.push_back(0);
    return;
  }
  digits.push_back(digits[digits.size() - 1]);
  for (size_t i = digits.size() - 2; i > 0; --i) {
    digits[i] = digits[i - 1];
  }
  digits[0] = 0;
}

BigInteger operator*(const BigInteger& left, const BigInteger& right) {
  BigInteger copy = left;
  copy *= right;
  return copy;
}

BigInteger& BigInteger::operator/=(const BigInteger& other) {
  if (other.digits.size() == 1 && other != 0) {
    long long divisor = other.digits[0];
    long long carry = 0;
    long long comp;
    for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
      comp = digits[i] + carry * kDigitLimit;
      digits[i] = comp / divisor;
      carry = comp % divisor;
    }
    removeLeadingZeros();
    sign = sign == other.sign ? Sign::PLUS : Sign::MINUS;
    return *this;
  }
  BigInteger copy_other = other >= 0 ? other : -other;
  BigInteger cur;
  for (int i = static_cast<int>(digits.size() - 1); i >= 0; --i) {
    cur.shift();
    cur.digits[0] = digits[i];
    cur.removeLeadingZeros();
    int x = 0, l = 0, r = kDigitLimit;
    while (l <= r) {
      int m = (l + r) / 2;
      BigInteger temp = copy_other * m;
      if (temp <= cur) {
        x = m;
        l = m + 1;
      } else
        r = m - 1;
    }
    digits[i] = x;
    cur -= copy_other * x;
  }
  removeLeadingZeros();
  sign = sign == other.sign ? Sign::PLUS : Sign::MINUS;
  normalizeSign();
  return *this;
}

BigInteger operator/(const BigInteger& left, const BigInteger& right) {
  BigInteger copy = left;
  copy /= right;
  return copy;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
  *this -= (*this / other) * other;
  return *this;
}

BigInteger operator%(const BigInteger& left, const BigInteger& right) {
  BigInteger copy = left;
  copy %= right;
  return copy;
}

BigInteger::operator bool() const {
  if (*this == 0) return false;
  return true;
}

BigInteger gcd(const BigInteger& left, const BigInteger& right) {
  BigInteger ans = 1;
  BigInteger copyLeft = left.abs();
  BigInteger copyRight = right.abs();

  while (copyLeft > 0 && copyRight > 0) {
    bool leftEven = copyLeft.even();
    bool rightEven = copyRight.even();
    if (leftEven && rightEven) {
      copyLeft /= 2;
      copyRight /= 2;
      ans *= 2;
    } else if (leftEven)
      copyLeft /= 2;
    else if (rightEven)
      copyRight /= 2;
    else
      copyRight < copyLeft ? copyLeft -= copyRight : copyRight -= copyLeft;
  }
  ans *= copyLeft != 0 ? copyLeft : copyRight;
  return ans;
}

class Rational {
  friend bool operator==(const Rational& left, const Rational& right);
  friend bool operator<(const Rational& left, const Rational& right);
  friend std::istream& operator>>(std::istream& in, Rational& rational);

 public:
  Rational();
  Rational(int number);
  Rational(const BigInteger& number);
  std::string toString() const;
  Rational operator-() const;
  Rational& operator+=(const Rational& other);
  Rational& operator++();
  Rational operator++(int);
  Rational& operator-=(const Rational& other);
  Rational& operator--();
  Rational operator--(int);
  Rational& operator*=(const Rational& other);
  Rational& operator/=(const Rational& other);
  std::string asDecimal(size_t precision) const;
  explicit operator double() const;

 private:
  BigInteger numerator;
  BigInteger denominator;
  void normSign();
  void reduce();
};

Rational::Rational() : numerator(), denominator() {}

Rational::Rational(int number) {
  numerator = number;
  denominator = 1;
}

Rational::Rational(const BigInteger& number) {
  numerator = number;
  denominator = 1;
}

void Rational::normSign() {
  if (denominator < 0) {
    numerator.changeSign();
    denominator.changeSign();
  }
}

void Rational::reduce() {
  BigInteger divider = gcd(numerator, denominator);
  if (divider != 1) {
    numerator /= divider;
    denominator /= divider;
  }
  normSign();
}

std::string Rational::toString() const {
  std::string str = numerator.toString();
  if (denominator != 1) str += '/' + denominator.toString();
  return str;
}

bool operator==(const Rational& left, const Rational& right) {
  if (left.numerator == right.numerator &&
      left.denominator == right.denominator)
    return true;
  if (left.numerator == 0 && right.numerator == 0) return true;
  return false;
}

bool operator<(const Rational& left, const Rational& right) {
  BigInteger l_num = left.numerator * right.denominator;
  BigInteger r_num = right.numerator * left.denominator;
  return l_num < r_num;
}

bool operator!=(const Rational& left, const Rational& right) {
  return !(left == right);
}

bool operator>(const Rational& left, const Rational& right) {
  return right < left;
}

bool operator<=(const Rational& left, const Rational& right) {
  return !(left > right);
}

bool operator>=(const Rational& left, const Rational& right) {
  return !(left < right);
}

Rational Rational::operator-() const {
  Rational copy;
  copy.numerator = -numerator;
  copy.denominator = denominator;
  return copy;
}

Rational& Rational::operator+=(const Rational& other) {
  BigInteger help = other.numerator * denominator;
  numerator *= other.denominator;
  denominator *= other.denominator;
  numerator += help;
  reduce();
  return *this;
}

Rational operator+(const Rational& left, const Rational& right) {
  Rational copy = left;
  copy += right;
  return copy;
}

Rational& Rational::operator++() {
  *this += 1;
  return *this;
}

Rational Rational::operator++(int) {
  Rational copy = *this;
  *this += 1;
  return copy;
}

Rational& Rational::operator-=(const Rational& other) {
  BigInteger help = other.numerator * denominator;
  numerator *= other.denominator;
  denominator *= other.denominator;
  numerator -= help;
  reduce();
  return *this;
}

Rational operator-(const Rational& left, const Rational& right) {
  Rational copy = left;
  copy -= right;
  return copy;
}

Rational& Rational::operator--() {
  *this -= 1;
  return *this;
}

Rational Rational::operator--(int) {
  Rational copy = *this;
  *this -= 1;
  return copy;
}

Rational& Rational::operator*=(const Rational& other) {
  if (other == 0) return *this = 0;
  if (other == 1) return *this;
  if (other == -1) {
    numerator.changeSign();
    return *this;
  }
  if (*this == 1) return *this = other;
  if (*this == -1) {
    *this = other;
    numerator.changeSign();
    return *this;
  }
  numerator *= other.numerator;
  denominator *= other.denominator;
  reduce();
  return *this;
}

Rational operator*(const Rational& left, const Rational& right) {
  Rational copy = left;
  copy *= right;
  return copy;
}

Rational& Rational::operator/=(const Rational& other) {
  numerator *= other.denominator;
  denominator *= other.numerator;
  reduce();
  return *this;
}

Rational operator/(const Rational& left, const Rational& right) {
  Rational copy = left;
  copy /= right;
  return copy;
}

std::string Rational::asDecimal(size_t precision) const {
  std::string decimal;
  std::string nulls;
  BigInteger help = (numerator >= 0) ? numerator : -numerator;
  for (size_t i = 0; i < precision; ++i) {
    help *= 10;
  }
  decimal += (help / denominator).toString();
  while (decimal.size() <= precision - nulls.size()) {
    nulls += '0';
  }
  decimal = ((numerator >= 0) ? nulls : '-' + nulls) + decimal;
  decimal.insert(decimal.size() - precision, ".");
  return decimal;
}

Rational::operator double() const { return std::stod(asDecimal(32)); }

std::ostream& operator<<(std::ostream& out, const Rational& rational) {
  out << rational.toString();
  return out;
}

std::istream& operator>>(std::istream& in, Rational& rational) {
  in >> rational.numerator;
  rational.denominator = 1;
  return in;
}
