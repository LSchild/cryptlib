//
// Created by Leonard on 7/23/26.
//

#include "backend/backend.h"
#include "backend/hexl_worker.h"

std::shared_ptr<MathWorker> SelectWorker(uint64_t modulus, uint64_t dimension) {

    // Future: select worker based on modulus & dimensio value

    return std::make_shared<HexlWorker>(modulus, dimension);
}