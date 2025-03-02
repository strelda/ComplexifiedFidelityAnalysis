#!/usr/bin/env julia

# Lipkin model parameters
const h_value        = BigFloat(1.0)    # strength of the interaction
const n_value        = 20                # Spin: j = n_value/2, matrix dimension = n_value + 1

# quench parameters
const lambda_initial = BigFloat(7.0)    # initial value of the quench parameter
const level_initial  = 1                # initial level of the quench
const lambda_final   = BigFloat(2.0)    # final value of the quench parameter

# time and beta ranges
const beta_range = (BigFloat("-2.5"), BigFloat("-1."))    # Real (β) range
const t_range    = (BigFloat("0.0"), BigFloat("20.0"))    # Imaginary (t) range

# precision settings
const precision = 64              # Digit precision settings: 128 above n_value=30, 256 for n_value=500
const recursion_steps = 2         # Maximum recursion depth for rectangle subdivision
const subd = 8                    # Number of subdivisions per edge

# Integration settings
const integration_order = 2                   # Used in quadgk
const integration_tolerance = 1e-2            # Used in quadgk
const winding_threshold = BigFloat("0.2")     # Lower threshold for nonzero winding

# resolution of the background image
const background_resolution = 50