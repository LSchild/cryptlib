//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_KEYS_H
#define LARGE_FUNCTIONS_KEYS_H

#include <string>
#include <cstdint>
#include <vector>
/**
 * Struct representing a key, any key
 */
struct GenericKey {

    GenericKey(const std::string& label, const std::vector<uint64_t>& key_vals);

    GenericKey(const std::string& label, const uint64_t* key_ptr, uint64_t n);

    GenericKey(const std::vector<uint64_t>& key_vals);

    GenericKey(const uint64_t* key_ptr, uint64_t n);

    [[nodiscard]] const uint64_t* GetKeyPtr() const {
        return m_key.data();
    }

    [[nodiscard]] const std::vector<uint64_t>& GetKey() const {
        return m_key;
    }

private:

    /* label associated with a key */
    std::string m_label;
    /* actual values of a key, flattened into a vector */
    std::vector<uint64_t> m_key;

};

#endif //LARGE_FUNCTIONS_KEYS_H
