"""Dispatch marker for a kinematics implementation."""
abstract type KinematicsMode end

"""Special-relativistic reaction kinematics."""
struct Relativistic <: KinematicsMode end
