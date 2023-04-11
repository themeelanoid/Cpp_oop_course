#include <cstring>
#include <iostream>

class String {
 public:
  String(char c);
  String(const char* c_string);
  String(size_t n, char c);
  String();
  ~String();
  String(const String& str);
  String& operator=(const String& str);
  char& operator[](size_t ind);
  const char& operator[](size_t ind) const;
  char& front();
  const char& front() const;
  char& back();
  const char& back() const;
  size_t length() const;
  void push_back(char symbol);
  void pop_back();
  String& operator+=(const String& str);
  size_t find(const String& substr) const;
  size_t rfind(const String& substr) const;
  String substr(size_t start, size_t count) const;
  bool empty() const;
  void clear();
  friend std::ostream& operator<<(std::ostream& out, const String& str);
  friend std::istream& operator>>(std::istream& in, String& str);

 private:
  char* string;
  size_t size;
  size_t volume;
  void exchange(String& str);
};

String::String(char c) : string(new char[1]), size(1), volume(1) {
  string[0] = c;
}

void String::exchange(String& str) {
  std::swap(str.size, size);
  std::swap(str.volume, volume);
  std::swap(str.string, string);
}

String::String(const char* c_string)
    : string(new char[strlen(c_string) * 2]),
      size(strlen(c_string)),
      volume(strlen(c_string) * 2) {
  for (size_t i = 0; i < size; ++i) {
    string[i] = c_string[i];
  }
}

String::String(size_t n, char c)
    : string(new char[n * 2]), size(n), volume(n * 2) {
  memset(string, c, n);
}

String::String() : string(new char), size(0), volume(0) {}

String::~String() { delete[] string; }

String::String(const String& str)
    : string(new char[str.volume]), size(str.size), volume(str.volume) {
  memcpy(string, str.string, sizeof(char) * (str.size + 1));
}

String& String::operator=(const String& str) {
  String copy = str;
  exchange(copy);
  return *this;
}

char& String::operator[](size_t ind) { return string[ind]; }

const char& String::operator[](size_t ind) const { return string[ind]; }

char& String::front() { return string[0]; }

const char& String::front() const { return string[0]; }

char& String::back() { return string[size - 1]; }

const char& String::back() const { return string[size - 1]; }

size_t String::length() const { return size; }

void String::push_back(char symbol) {
  if (size == volume - 1) {
    char* newstring = new char[volume * 2];
    memcpy(newstring, string, sizeof(char) * (size + 1));
    delete[] string;
    string = newstring;
    volume *= 2;
  }
  string[size] = symbol;
  ++size;
}

void String::pop_back() {
  if (size > 0) --size;
}

String& String::operator+=(const String& str) {
  if (volume <= size + str.size) {
    char* newstring = new char[(size + str.size) * 2];
    memcpy(newstring, string, sizeof(char) * size);
    delete[] string;
    string = newstring;
    volume = (size + str.size) * 2;
  }
  for (size_t i = 0; i < str.size; ++i) {
    string[size + i] = str[i];
  }
  size += str.size;
  return *this;
}

String operator+(const String& str1, const String& str2) {
  String newstring(str1);
  newstring += str2;
  return newstring;
}

size_t String::find(const String& substr) const {
  for (size_t i = 0; i < size - substr.size; ++i) {
    size_t j = 0;
    while (j < substr.size && string[i + j] == substr[j]) {
      if (j == substr.size - 1) return i;
      ++j;
    }
  }
  return this->length();
}

size_t String::rfind(const String& substr) const {
  for (size_t i = size - 1; i >= substr.size; --i) {
    size_t j = substr.size;
    while (j >= 1 && string[i - substr.size + j] == substr[j - 1]) {
      if (j == 1) return i - substr.size + 1;
      --j;
    }
  }
  return this->length();
}

String String::substr(size_t start, size_t count) const {
  String substr(count, ' ');
  memcpy(substr.string, string + start, sizeof(char) * size);
  return substr;
}

bool String::empty() const { return size == 0; }

void String::clear() { size = 0; }

std::ostream& operator<<(std::ostream& out, const String& str) {
  for (size_t i = 0; i < str.size; ++i) {
    out << str.string[i];
  }
  return out;
}

std::istream& operator>>(std::istream& in, String& str) {
  char c;
  while (isspace(in.peek())) {
    in.get();
  }
  while (!isspace(in.peek()) && in >> c) {
    str.push_back(c);
  }
  return in;
}

bool operator==(const String& left, const String& right) {
  if (left.length() != right.length()) return false;
  for (size_t i = 0; i < left.length(); ++i) {
    if (left[i] != right[i]) return false;
  }
  return true;
}