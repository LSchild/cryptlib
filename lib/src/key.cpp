//
// Created by leonard on 5/27/26.
//

#include "keys.h"

GenericKey::GenericKey(const std::string &label, const std::vector<uint64_t>& key_vals) : m_label(label) {
    m_key.resize(key_vals.size());
    std::copy(key_vals.begin(), key_vals.end(), m_key.begin());
}

GenericKey::GenericKey(const std::string &label, const uint64_t *key_ptr, uint64_t n) : m_label(label) {
    m_key.resize(n);
    std::copy(key_ptr, key_ptr + n, m_key.data());
}

GenericKey::GenericKey(const uint64_t *key_ptr, uint64_t n) : m_label("EMPTY_FROM_PTR") {
    m_key.resize(n);
    std::copy(key_ptr, key_ptr + n, m_key.data());
}

GenericKey::GenericKey(const std::vector<uint64_t>& key_vals) : m_label("EMPTY_FROM_VEC") {
    m_key.resize(key_vals.size());
    std::copy(key_vals.begin(), key_vals.end(), m_key.begin());
}