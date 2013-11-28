disp('Position errors')
[results.M1withoutVel.kita_rmse.pos]
[results.M3withoutVel.kita_rmse.pos]

disp('Velocity errors')
[results.M1withoutVel.kita_rmse.vel]
[results.M3withoutVel.kita_rmse.vel]

disp('Mean errors: intrinsic model, unobserved velocity')
mean([results.M1withoutVel.filt_rmse.pos])
mean([results.M1withoutVel.filt_rmse.vel])
mean([results.M1withoutVel.kita_rmse.pos])
mean([results.M1withoutVel.kita_rmse.vel])

disp('Mean errors: Cartesian model, unobserved velocity')
mean([results.M3withoutVel.filt_rmse.pos])
mean([results.M3withoutVel.filt_rmse.vel])
mean([results.M3withoutVel.kita_rmse.pos])
mean([results.M3withoutVel.kita_rmse.vel])

disp('Error standard deviation: intrinsic model, unobserved velocity')
std([results.M1withoutVel.filt_rmse.pos])
std([results.M1withoutVel.filt_rmse.vel])
std([results.M1withoutVel.kita_rmse.pos])
std([results.M1withoutVel.kita_rmse.vel])

disp('Error standard deviation: Cartesian model, unobserved velocity')
std([results.M3withoutVel.filt_rmse.pos])
std([results.M3withoutVel.filt_rmse.vel])
std([results.M3withoutVel.kita_rmse.pos])
std([results.M3withoutVel.kita_rmse.vel])