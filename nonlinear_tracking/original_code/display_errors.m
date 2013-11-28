disp('Position errors')
[results.M1withoutVel.kita_rmse.pos]
[results.M1withVel.kita_rmse.pos]
[results.M2withoutVel.kita_rmse.pos]
[results.M2withVel.kita_rmse.pos]

disp('Velocity errors')
[results.M1withoutVel.kita_rmse.vel]
[results.M1withVel.kita_rmse.vel]
[results.M2withoutVel.kita_rmse.vel]
[results.M2withVel.kita_rmse.vel]

disp('Mean errors: basic model, unobserved velocity')
mean([results.M1withoutVel.filt_rmse.pos])
mean([results.M1withoutVel.filt_rmse.vel])
mean([results.M1withoutVel.kita_rmse.pos])
mean([results.M1withoutVel.kita_rmse.vel])

disp('Mean errors: basic model, observed velocity')
mean([results.M1withVel.filt_rmse.pos])
mean([results.M1withVel.filt_rmse.vel])
mean([results.M1withVel.kita_rmse.pos])
mean([results.M1withVel.kita_rmse.vel])

disp('Mean errors: augmented model, unobserved velocity')
mean([results.M2withoutVel.filt_rmse.pos])
mean([results.M2withoutVel.filt_rmse.vel])
mean([results.M2withoutVel.kita_rmse.pos])
mean([results.M2withoutVel.kita_rmse.vel])

disp('Mean errors: augmented model, observed velocity')
mean([results.M2withVel.filt_rmse.pos])
mean([results.M2withVel.filt_rmse.vel])
mean([results.M2withVel.kita_rmse.pos])
mean([results.M2withVel.kita_rmse.vel])

disp('Error standard deviation: basic model, unobserved velocity')
std([results.M1withoutVel.filt_rmse.pos])
std([results.M1withoutVel.filt_rmse.vel])
std([results.M1withoutVel.kita_rmse.pos])
std([results.M1withoutVel.kita_rmse.vel])

disp('Error standard deviation: basic model, observed velocity')
std([results.M1withVel.filt_rmse.pos])
std([results.M1withVel.filt_rmse.vel])
std([results.M1withVel.kita_rmse.pos])
std([results.M1withVel.kita_rmse.vel])

disp('Error standard deviation: augmented model, unobserved velocity')
std([results.M2withoutVel.filt_rmse.pos])
std([results.M2withoutVel.filt_rmse.vel])
std([results.M2withoutVel.kita_rmse.pos])
std([results.M2withoutVel.kita_rmse.vel])

disp('Error standard deviation: augmented model, observed velocity')
std([results.M2withVel.filt_rmse.pos])
std([results.M2withVel.filt_rmse.vel])
std([results.M2withVel.kita_rmse.pos])
std([results.M2withVel.kita_rmse.vel])