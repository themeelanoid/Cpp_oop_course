#include <iostream>

template <size_t N>
class StackStorage {
 public:
  StackStorage() = default;
  uint8_t* allocate(size_t n, const size_t align);
  StackStorage(const StackStorage& other) = delete;
  void deallocate(uint8_t* ptr, size_t n, const size_t align);
  StackStorage& operator=(const StackStorage& other) = delete;

 private:
  uint8_t arr[N];
  uint8_t* pointer = arr;
};

template <size_t N>
uint8_t* StackStorage<N>::allocate(size_t n, const size_t align) {
  pointer += n + (n - reinterpret_cast<size_t>(pointer)) % align;
  return pointer - n;
}

template <size_t N>
void StackStorage<N>::deallocate(uint8_t* ptr, size_t n, const size_t align) {
  n += align - ((reinterpret_cast<size_t>(pointer) + n) % align);
  if (ptr + n == pointer) {
    pointer = ptr;
  }
}

template <typename T, size_t N>
class StackAllocator {
 public:
  typedef T value_type;
  StackAllocator() = default;
  ~StackAllocator() = default;
  StackAllocator(StackStorage<N>& storage) : storage(&storage) {}
  template <typename U>
  struct rebind {
    typedef StackAllocator<U, N> other;
  };
  T* allocate(size_t n) {
    return reinterpret_cast<T*>(storage->allocate(n * sizeof(T), alignof(T)));
  }
  void deallocate(T* p, size_t n) {
    storage->deallocate(reinterpret_cast<uint8_t*>(p), n * sizeof(T),
                        alignof(T));
  }
  template <typename U>
  StackAllocator(const StackAllocator<U, N>& other) : storage(other.storage) {}
  StackStorage<N>* storage;
};

template <typename T, typename Alloc = std::allocator<T>>
class List {
 public:
  template <bool isConst>
  class Iter;

  typedef Iter<false> iterator;
  iterator begin() { return iterator(fakeNode.next); }
  iterator end() { return iterator(&fakeNode); }

  typedef Iter<true> const_iterator;
  const_iterator begin() const { return const_iterator(fakeNode.next); }
  const_iterator end() const { return const_iterator(&fakeNode); }
  const_iterator cbegin() const { return const_iterator(fakeNode.next); }
  const_iterator cend() const { return const_iterator(&fakeNode); }

  typedef std::reverse_iterator<Iter<false>> reverse_iterator;
  reverse_iterator rbegin() { return std::reverse_iterator(end()); }
  reverse_iterator rend() { return std::reverse_iterator(begin()); }

  typedef std::reverse_iterator<Iter<true>> const_reverse_iterator;
  const_reverse_iterator rend() const { return std::reverse_iterator(begin()); }
  const_reverse_iterator crend() const {
    return std::reverse_iterator(begin());
  }
  const_reverse_iterator rbegin() const {
    return std::reverse_iterator(cend());
  }
  const_reverse_iterator crbegin() const {
    return std::reverse_iterator(cend());
  }

  List(int n, const T& value, const Alloc& alloc);
  explicit List(int n);
  explicit List(Alloc& alloc) : alloc(alloc) {}
  List(int n, const T& value);
  List(int n, Alloc& alloc);
  Alloc get_allocator() const { return alloc; }
  List() = default;
  ~List() { clear(); }
  size_t size() const { return sz; };
  List(const List& other);
  List& operator=(const List& other);

  void push_back(const T& value);
  void push_front(const T& value);
  void erase(const_iterator iter);
  void pop_back();
  void insert(const_iterator iter, const T& value);
  void pop_front();
  void clear();

 private:
  struct BaseNode {
    BaseNode* prev;
    BaseNode* next;
    BaseNode() : prev(this), next(this) {}
    BaseNode(BaseNode* p, BaseNode* n) : prev(p), next(n) {}
  };

  struct Node : BaseNode {
    Node() = default;
    T value;
    Node(const T& value, BaseNode* prev, BaseNode* next)
        : BaseNode(prev, next), value(value) {}
  };

  void connectNodes(BaseNode* node, BaseNode* p, BaseNode* n) {
    node->next->prev = p;
    node->prev->next = n;
  }

  typedef std::allocator_traits<
      typename std::allocator_traits<Alloc>::template rebind_alloc<Node>>
      NodeAllocTraits;

  typename std::allocator_traits<Alloc>::template rebind_alloc<Node> alloc;
  size_t sz = 0;
  mutable BaseNode fakeNode;
};

template <typename T, typename Alloc>
template <bool isConst>
class List<T, Alloc>::Iter {
 public:
  BaseNode* iter;

  Iter() = default;
  Iter(BaseNode* node) : iter(node){};
  Iter(const Iter<false>& other) : iter(other.iter){};

  typedef std::bidirectional_iterator_tag iterator_category;
  typedef T value_type;
  typedef T* pointer;
  typedef std::ptrdiff_t difference_type;
  typedef std::conditional_t<isConst, const T&, T&> reference;

  Iter& operator++() {
    iter = iter->next;
    return *this;
  }
  bool operator==(const Iter<isConst>& other) const {
    return (iter == other.iter);
  }
  Iter& operator--() {
    iter = iter->prev;
    return *this;
  }
  std::conditional_t<isConst, const T&, T&> operator*() const {
    return reinterpret_cast<Node*>(iter)->value;
  }
  bool operator!=(const Iter<isConst>& other) const {
    return !(*this == other);
  }
  Iter operator++(int) {
    Iter copy = *this;
    iter = iter->next;
    return copy;
  }
  Iter operator--(int) {
    Iter copy = *this;
    iter = iter->prev;
    return copy;
  }
  std::conditional_t<isConst, const T*, T*> operator->() const {
    return &(reinterpret_cast<Node*>(iter)->value);
  }
};

template <typename T, typename Alloc>
void List<T, Alloc>::insert(List::const_iterator iter, const T& value) {
  Node* node = NodeAllocTraits::allocate(alloc, 1);
  try {
    NodeAllocTraits::construct(alloc, node, value, iter.iter->prev, iter.iter);
    ++sz;
    connectNodes(node, node, node);
  } catch (...) {
    NodeAllocTraits::deallocate(alloc, node, 1);
    throw;
  }
}

template <typename T, typename Alloc>
void List<T, Alloc>::erase(List::const_iterator iter) {
  connectNodes(iter.iter, iter.iter->prev, iter.iter->next);
  --sz;
  NodeAllocTraits::destroy(alloc, static_cast<Node*>(iter.iter));
  NodeAllocTraits::deallocate(alloc, static_cast<Node*>(iter.iter), 1);
}

template <typename T, typename Alloc>
void List<T, Alloc>::pop_back() {
  erase(std::prev(end()));
}

template <typename T, typename Alloc>
void List<T, Alloc>::clear() {
  while (sz > 0) {
    pop_back();
  }
}

template <typename T, typename Alloc>
void List<T, Alloc>::push_back(const T& value) {
  insert(end(), value);
}

template <typename T, typename Alloc>
List<T, Alloc>::List(int n, const T& value) {
  try {
    while (n--) {
      push_front(value);
    }
  } catch (...) {
    clear();
  }
}

template <typename T, typename Alloc>
List<T, Alloc>::List(int n) {
  try {
    while (n--) {
      Node* node = NodeAllocTraits::allocate(alloc, 1);
      try {
        NodeAllocTraits::construct(alloc, node);
        node->next = end().iter;
        node->prev = end().iter->prev;
      } catch (...) {
        NodeAllocTraits::deallocate(alloc, node, 1);
        throw;
      }
      ++sz;
      connectNodes(node, node, node);
    }
  } catch (...) {
    clear();
  }
}

template <typename T, typename Alloc>
void List<T, Alloc>::pop_front() {
  erase(begin());
}

template <typename T, typename Alloc>
List<T, Alloc>::List(const List& other)
    : alloc(NodeAllocTraits::select_on_container_copy_construction(
          other.get_allocator())) {
  try {
    for (auto it = other.begin(); it != other.end(); ++it) {
      push_back(*it);
    }
  } catch (...) {
    clear();
    throw;
  }
}

template <typename T, typename Alloc>
void List<T, Alloc>::push_front(const T& value) {
  insert(begin(), value);
}

template <typename T, typename Alloc>
List<T, Alloc>::List(int n, const T& value, const Alloc& alloc) : alloc(alloc) {
  try {
    while (n--) {
      push_back(value);
    }
  } catch (...) {
    clear();
  }
}

template <typename T, typename Alloc>
List<T, Alloc>::List(int n, Alloc& allocator) : alloc(allocator) {
  try {
    while (n--) {
      Node* node = NodeAllocTraits::allocate(alloc, 1);
      try {
        NodeAllocTraits::construct(alloc, node);
        node->next = end().iter;
        node->prev = end().iter->prev;
      } catch (...) {
        NodeAllocTraits::deallocate(alloc, node, 1);
        throw;
      }
      ++sz;
      connectNodes(node, node, node);
    }
  } catch (...) {
    clear();
  }
}

template <typename T, typename Alloc>
List<T, Alloc>& List<T, Alloc>::operator=(const List& other) {
  if (this != &other) {
    int prev_size = static_cast<int>(sz);
    int added = 0;
    try {
      for (auto it = other.begin(); it != other.end(); ++it) {
        push_back(*it);
        ++added;
      }
    } catch (...) {
      while (added--) {
        pop_back();
      }
      throw;
    }
    if (NodeAllocTraits::propagate_on_container_copy_assignment::value) {
      alloc = other.alloc;
    }
    while (prev_size--) {
      pop_front();
    }
  }
  return *this;
}