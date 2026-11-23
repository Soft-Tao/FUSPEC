const RNG_GOLDEN_GAMMA = UInt64(0x9e3779b97f4a7c15)
const RNG_SAMPLE_MIX = UInt64(0xd2b74407b1ce6e93)
const RNG_DIMENSION_MIX = UInt64(0xca5a826395121157)

@inline function mix64(value::UInt64)
    value = xor(value, value >> 30)
    value *= UInt64(0xbf58476d1ce4e5b9)
    value = xor(value, value >> 27)
    value *= UInt64(0x94d049bb133111eb)
    return xor(value, value >> 31)
end

"""Counter-based random bits, independent of call and thread order."""
@inline function counter_u64(seed::UInt64, sample_index::UInt64,
                             dimension::UInt32)
    counter = seed + RNG_GOLDEN_GAMMA
    counter = xor(counter, (sample_index + UInt64(1)) * RNG_SAMPLE_MIX)
    counter = xor(counter, (UInt64(dimension) + UInt64(1)) *
                           RNG_DIMENSION_MIX)
    return mix64(counter)
end

counter_u64(seed::Integer, sample_index::Integer, dimension::Integer) =
    counter_u64(UInt64(seed), UInt64(sample_index), UInt32(dimension))

"""Uniform random value in [0, 1) for a fixed sample and random dimension."""
@inline counter_uniform(seed::UInt64, sample_index::UInt64,
                        dimension::UInt32) =
    Float64(counter_u64(seed, sample_index, dimension) >> 11) * 0x1.0p-53

counter_uniform(seed::Integer, sample_index::Integer, dimension::Integer) =
    counter_uniform(UInt64(seed), UInt64(sample_index), UInt32(dimension))
