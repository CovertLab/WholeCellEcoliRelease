This is a minimal example of a metabolic model that is driven by the difference between target (homeostatic) and actual metabolite concentrations, using the formalism of a proportional-integral (PI) controller.  It is *not* a complete metabolic model, but rather a model of the behavior of a metabolic model, which in the context of a whole-cell model would need to be incorporated with a dynamic flux-balance analysis (dFBA)-like solver to determine what metabolism actually accomplishes (and how).

# Background

The motivation is subtle.  Currently, in the *E. coli* whole-cell model, metabolism is (under sufficient conditions) designed to solve

```
0 = c_final - c_target
```

for each metabolite, where `c_target` is some target concentration, and `c_final` is the 'partitioned concentration' after metabolism has run.  That is, `c_final = c_partitioned + v * dt` where `c_partitioned` is, of course, what is partitioned, `dt` is the length of a time-step, and `v` is the 'velocity' or flux (rate) at which that molecule is produced.  Rearranging, we obtain

```
v = (c_target - c_partitioned) / dt
```

Now, for reasons too complicated to go into here, removing this 'partitioning' idea is desirable.  Naively we replace `c_partitioned` - which is artifical and process-specific - with `c_initial`, the real, global concentration at the beginning of the time step:

```
v = (c_target - c_initial) / dt
```

This sounds fine, and in practice it sort of works, but due to the discretization of time, the target is never actually hit (we always fall short).  This issue gets worse with longer time-steps, to a point where the system becomes numerically unstable.

This is resolved by observing that what we have built is effectively a proportional controller, i.e. we are setting `v` as proportional to our error (target less initial).  The proportionality constant is the inverse of our step size, but in reality we could choose *any* proportionality, with the standard caveats; too large, and we will catastrophically overshoot our target; too small, and we will never reach our target.

One solution - which works quite well - is to set `v` based not only on the instantaneous error (which may be too short-sighted), but also on the *accumulated* error over all time steps.  This is known as integral control (since we are *integrating* the error over all time); it has its own issues, but with proper tuning it can address the issues that crop up with proportional control.  Roughly, the algorithm is

```
accumulated_error = 0

while simulating:
	error = c_target - c_initial

	v = dt * (k_prop * error + k_integral * accumulated_error)

	accumulated_error += error
```

This is exactly what I implement in this demonstration.  My context is a two-metabolite system of ATP and one amino acid (alanine), both of which are utilized by a simplified protein translation model.  Through a series of interventions I go from the naive model to the updated PI metabolic model.  Notably this enables longer time-steps, and when combined with a slightly more complex translation model, also allows for realistic metabolic down-shifts.

# Demonstrations

## Demonstration 1

The naive implementation, using a step-size of 1.0 seconds (equivalent to our current model).  Note that molecule abundances are always short of their targets (dashed lines).

## Demonstration 2

The same as #1, but with a time-step length of 10.0 seconds, illustrating the basic pathology associated with partitioning.

## Demonstration 3

Back to a time-step length of 1.0 seconds, but now with integral control.  Note the short acclimation period where the integral control error is accumulating, after which the target abundances are hit more or less exactly.

## Demonstration 4

The same as #3, but with the accumulated error 'bootstrapped' (using values from a prior simulation) to start near the right target value, largely eliminating the acclimation period.  This bootstrap value could probably be determined analytically - I haven't tried to do so.

## Demonstration 5

The same as #4, except that I have increased the time-step length to 10.0 seconds.  Contrast with #2.  (Note that bootstrapping is needed to avoid numerical instability - in practice, slowly ramping up the time-step length during the acclimation period would also work.)

## Demonstration 6

The same as #5, except that the rate of alanine production is artificially constrained to about 90% of its anticipated demand (to simulate a metabolic downshift).  Note the linear fall-off until zero abundance is reached, after which time a pathological behavior emerges.  This pathology is a consequence of the fact that the activity of the translation model does not slow down until it hits the hard limit of zero available molecules.

## Demonstration 7

The same as #6, except that translation is now modeled such that it begins to slow down once either of its substrates fall to about 1% of their target concentrations.  Note the smooth fall-off to a new equilibrium value for both the amino acid abundance as well as the overall translation rate.

# Biological notes

The proportional controller may be rationalized as metabolic (e.g. allosteric) regulation.  The integral controller can be rationalized as a build up of potential (e.g. metabolite pooling) in the network, also consequent of regulation.

The last element introduced (the rate-limiting terms in the translation model) can be thought of as kinetic limitations on tRNA charging.

# Speculative WCM implementation

Bringing this to the whole-cell model would require redesigning metabolism to accept target fluxes rather than target concentrations - not a big change.

It would also require the introduction of artificial state variables (the accumulated error for the integral controller), which is unfortunate but I believe is well justified (see last section).  Notably these are floating-point variables, not non-negative integers, so this couldn't just be attached to `BulkMolecules`.  A new state object for this sort of variable would be ideal - in the intervening time, a `Listener` could be used.

Finally, the proportional and integral control constants would need to be tuned.  Apart from hand-tuning, there is extensive literature describing strategies for tuning these constants.  I would first start with the values I selected and see just how bad those are.

The rate-limiting terms in the kinetic model aren't eminently needed but they would make the model output more reasonable.  The ongoing work on kinetic tRNA charging will eventually address this.
