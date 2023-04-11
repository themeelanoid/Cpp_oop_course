#include <vector>

template <typename T>
class Deque {
 private:
  std::vector<T*> pointers;
  static const size_t block_size = 64;
  size_t start_block_number = 0;
  size_t start_index = 0;
  size_t end_block_number = 0;
  size_t end_index = 0;
  size_t sz = 0;

 public:
  void swap(Deque& deque) {
    std::swap(start_block_number, deque.start_block_number);
    std::swap(start_index, deque.start_index);
    std::swap(end_block_number, deque.end_block_number);
    std::swap(end_index, deque.end_index);
    std::swap(sz, deque.sz);
    std::swap(pointers, deque.pointers);
  }

  Deque() {
    pointers.resize(2);
    size_t block_number = pointers.size() / 2;
    start_block_number = block_number;
    start_index = block_size / 2;
    end_block_number = block_number;
    end_index = block_size / 2 - 1;
    pointers[0] = reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]);
    pointers[1] = reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]);
  }

  explicit Deque(const int& n, const T& value = T()) {
    size_t blocksCount = static_cast<size_t>(n) / block_size + 1;
    start_block_number = 0;
    start_index = 0;
    end_block_number = blocksCount - 1;
    end_index = (n - 1) % block_size;
    sz = n;
    pointers.assign(blocksCount,
                    reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]));
    for (size_t i = 0; i < blocksCount; ++i) {
      for (size_t j = 0; j < block_size; ++j) {
        if (i < blocksCount - 1 || j <= end_index) {
          new (pointers[i] + j) T(value);
        }
      }
    }
  }

  Deque(const Deque& deque) : Deque() {
    pointers.assign(deque.pointers.size(),
                    reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]));
    start_block_number = deque.start_block_number;
    end_block_number = deque.end_block_number;
    for (size_t i = 0; i < pointers.size(); ++i) {
      if (i >= start_block_number && i <= end_block_number) {
        for (size_t j = 0; j < block_size; ++j) {
          new (pointers[i] + j) T(deque.pointers[i][j]);
        }
      }
    }
    start_index = deque.start_index;
    end_index = deque.end_index;
    sz = deque.sz;
  }

  Deque& operator=(const Deque& deque) {
    Deque copy = Deque(deque);
    swap(copy);
    return *this;
  }

  size_t size() const { return sz; }

  T& operator[](size_t i) {
    size_t position = start_block_number * block_size + start_index + i;
    size_t block_number = position / block_size;
    size_t index = position % block_size;
    return pointers[block_number][index];
  }

  const T& operator[](size_t i) const {
    size_t position = start_block_number * block_size + start_index + i;
    size_t block_number = position / block_size;
    size_t index = position % block_size;
    return pointers[block_number][index];
  }

  T& at(size_t i) {
    if (i >= sz) throw std::out_of_range("Error: out of range");
    size_t position = start_block_number * block_size + start_index + i;
    size_t block_number = position / block_size;
    size_t index = position % block_size;
    return pointers[block_number][index];
  }

  const T& at(size_t i) const {
    if (i >= sz) throw std::out_of_range("Error: out of range");
    size_t position = start_block_number * block_size + start_index + i;
    size_t block_number = position / block_size;
    size_t index = position % block_size;
    return pointers[block_number][index];
  }

  void push_back(const T& value) {
    ++sz;
    if (end_block_number == pointers.size() - 1 &&
        end_index == block_size - 1) {
      pointers.resize(2 * pointers.size());
      for (size_t i = pointers.size() / 2; i < pointers.size(); ++i) {
        pointers[i] = reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]);
      }
    }

    if (end_index != block_size - 1) {
      ++end_index;
    } else {
      ++end_block_number;
      end_index = 0;
    }
    new (pointers[end_block_number] + end_index) T(value);
  }

  void pop_back() {
    --sz;
    (pointers[end_block_number] + end_index)->~T();
    if (end_index != 0) {
      --end_index;
    } else {
      --end_block_number;
      end_index = block_size - 1;
    }
  }

  void push_front(const T& value) {
    ++sz;
    if (start_block_number == 0 && start_index == 0) {
      pointers.resize(2 * pointers.size());
      start_block_number = pointers.size() / 2;
      end_block_number += pointers.size() / 2;
      for (size_t i = 0; i < pointers.size() / 2; ++i) {
        pointers[i + pointers.size() / 2] = pointers[i];
        pointers[i] = reinterpret_cast<T*>(new int8_t[block_size * sizeof(T)]);
      }
    }
    if (start_index != 0) {
      --start_index;
    } else {
      --start_block_number;
      start_index = block_size - 1;
    }
    new (pointers[start_block_number] + start_index) T(value);
  }

  void pop_front() {
    --sz;
    (pointers[start_block_number] + start_index)->~T();
    if (start_index != block_size - 1) {
      ++start_index;
    } else {
      ++start_block_number;
      start_index = 0;
    }
  }

  template <bool is_const>
  struct common_iterator {
    typename std::vector<T*>::const_iterator pointer;
    size_t index;

    common_iterator(typename std::vector<T*>::const_iterator ptr, size_t ind) {
      pointer = ptr;
      index = ind;
    }

    int operator-(const common_iterator<is_const>& other) const {
      int ans = (pointer - other.pointer) * block_size;
      ans += static_cast<int>(index) - static_cast<int>(other.index);
      return ans;
    }

    bool operator<(const common_iterator<is_const>& other) const {
      return (*this - other) < 0;
    }

    common_iterator<is_const> operator+(const int& n) const {
      return common_iterator<is_const>(pointer + (index + n) / block_size,
                                       (index + n) % block_size);
    }

    bool operator==(const common_iterator<is_const>& other) const {
      return *this - other == 0;
    }

    common_iterator& operator++() {
      if (index < block_size - 1) {
        ++index;
      } else {
        ++pointer;
        index = 0;
      }
      return *this;
    }

    bool operator<=(const common_iterator<is_const>& other) const {
      return *this < other || *this == other;
    }

    bool operator>=(const common_iterator<is_const>& other) const {
      return other <= *this;
    }

    common_iterator<is_const> operator-(const int& n) const {
      auto new_n = static_cast<size_t>(n);
      if (index >= new_n) {
        return common_iterator<is_const>(pointer + (index - new_n) / block_size,
                                         (index - new_n) % block_size);
      }
      int minus_block = (index - new_n) / block_size - 1;
      int new_elem_index = (index - new_n) % block_size + block_size;
      return common_iterator<is_const>(pointer + minus_block, new_elem_index);
    }

    common_iterator& operator--() {
      if (index > 0) {
        --index;
      } else {
        index = block_size - 1;
        --pointer;
      }
      return *this;
    }

    bool operator!=(const common_iterator<is_const>& other) const {
      return !(*this == other);
    }

    bool operator>(const common_iterator<is_const>& other) const {
      return other < *this;
    }

    std::conditional_t<is_const, const T&, T&> operator*() const {
      auto tmp = *pointer + index;
      return *tmp;
    }

    std::conditional_t<is_const, const T*, T*> operator->() const {
      auto tmp = *pointer + index;
      return tmp;
    }
  };

  using iterator = common_iterator<false>;
  using const_iterator = const common_iterator<true>;

  void insert(iterator it, const T& value) {
    push_back(value);
    for (auto i = --end(); i != it; --i) {
      std::iter_swap(i, i - 1);
    }
  }

  iterator begin() {
    return iterator(pointers.begin() + start_block_number, start_index);
  }

  iterator end() {
    auto copy = iterator(pointers.begin() + end_block_number, end_index);
    ++copy;
    return copy;
  }

  const_iterator begin() const {
    return const_iterator(pointers.begin() + start_block_number, start_index);
  }

  const_iterator end() const {
    auto copy = const_iterator(pointers.begin() + end_block_number, end_index);
    return ++copy;
  }

  const_iterator cbegin() const {
    return const_iterator(pointers.begin() + start_block_number, start_index);
  }

  const_iterator cend() const {
    auto copy = const_iterator(pointers.begin() + end_block_number, end_index);
    return ++copy;
  }

  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  std::reverse_iterator<iterator> rbegin() {
    return std::reverse_iterator(end());
  }

  std::reverse_iterator<iterator> rend() {
    return std::reverse_iterator(begin());
  }

  std::reverse_iterator<const_iterator> rbegin() const {
    return std::reverse_iterator(cend());
  }

  std::reverse_iterator<const_iterator> rend() const {
    return std::reverse_iterator(cbegin());
  }

  std::reverse_iterator<iterator> crbegin() const {
    return std::reverse_iterator(cend());
  }

  std::reverse_iterator<const_iterator> crend() const {
    return std::reverse_iterator(cbegin());
  }

  void erase(iterator it) {
    for (auto i = it; i != ++begin(); --i) {
      std::iter_swap(i - 1, i);
    }
    pop_front();
  }
};
