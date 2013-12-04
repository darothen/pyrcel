module shr_kind_mod
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer
   integer,parameter :: SHR_KIND_CS = 80                     ! short char
   integer,parameter :: SHR_KIND_CL = 256                    ! long char
   integer,parameter :: SHR_KIND_CX = 512                    ! extra-long char

end module shr_kind_mod

module MITaer_activation

use shr_kind_mod, only: r8=>shr_kind_r8

implicit none
public
save

contains

!=======================================================================
subroutine modal_smax( N, mu, sigma, kappa, V, T, P, smax )
!=======================================================================
! Purpose:
!  Computes the supersaturaton max achieved in a given updraft with a
!  prescribed aerosol.
!
! Developed using a polynomial chaos fit using quadrature based on a 
! 3rd-order Smolyak sparse grid (non-nested). All parameters were model-
! ed as uniform random variables with equal likelihood across an 
! estimate of their likely range. An exception is for the number 
! concentration, which spans several orders of magnitude and was modeled
! as a very wide lognormal distribution.
!
! If a negative value is returned, then either the given parameters 
! were outside of the fitted range or there was an error.
!
! Author: Daniel Rothenberg, 12/02/2013
!-----------------------------------------------------------------------

   implicit NONE

!---------------------------- Arguments --------------------------------
   integer :: nmode ! number of modes

   real(r8), intent(in) :: N     ! num con, cm^-3
   real(r8), intent(in) :: mu    ! mean r, micron 
   real(r8), intent(in) :: sigma ! shape
   real(r8), intent(in) :: kappa ! hygroscopicity
   real(r8), intent(in) :: V     ! Updraft speed, m/s
   real(r8), intent(in) :: T     ! Temperature, K
   real(r8), intent(in) :: P     ! Pressure, Pa

   real(r8), intent(out) :: smax ! maximum supersaturation
!------------------------ Local work space  ----------------------------
   real(r8), parameter :: N_lo = 1.0_r8, N_hi = 1.0e6_r8
   real(r8), parameter :: mu_lo = 0.01_r8, mu_hi = 0.25_r8
   real(r8), parameter :: sigma_lo = 1.2_r8, sigma_hi = 3.0_r8
   real(r8), parameter :: kappa_lo = 0.1_r8, kappa_hi = 1.2_r8
   real(r8), parameter :: V_lo = 0.0_r8, V_hi = 4.0_r8
   real(r8), parameter :: T_lo = 235.0_r8, T_hi = 310.0_r8
   real(r8), parameter :: P_lo = 50000.0_r8, P_hi = 105000.0_r8

   real(r8), parameter :: N_lambda = 7.5, N_zeta = 1.95

   real(r8), parameter :: rmin = 1.e-30_r8

   real(r8) :: N_map, mu_map, sigma_map, kappa_map, V_map, T_map, P_map
   logical :: out_of_range
!-----------------------------------------------------------------------

   ! Check that the parameters are okay and map them into standard
   ! random variables
   write (*, *) N, mu, sigma, kappa, V, T, P
   smax = -1.0_r8
   out_of_range = .false.
   if ((( mu < mu_lo ) .or. ( mu > mu_hi )) .or. &
       (( N < N_lo ) .or. ( N > N_hi )) .or. &
       (( sigma < sigma_lo ) .or. ( sigma > sigma_hi )) .or. &
       (( kappa < kappa_lo ) .or. ( kappa > kappa_hi )) .or. &
       (( V < V_lo ) .or. ( V > V_hi )) .or. &
       (( P < P_lo ) .or. ( P > P_hi )) .or. &
       (( T < T_lo ) .or. ( T > T_hi ))) return

   ! Map the parameter into the appropriate space
   N_map = lognorm_to_norm(N, N_lambda, N_zeta)

   mu_map = uni_to_norm(mu, mu_lo, mu_hi)
   kappa_map = uni_to_norm(kappa, kappa_lo, kappa_hi)
   sigma_map = uni_to_norm(sigma, sigma_lo, sigma_hi)
   V_map = uni_to_norm(V, V_lo, V_hi)
   T_map = uni_to_norm(T, T_lo, T_hi)
   P_map = uni_to_norm(P, P_lo, P_hi)

   smax = 6.4584111537e-6*N_map**6 - &
         2.5994976288e-5*N_map**5 - &
         1.7065251097e-7*N_map**4*P_map**2 + &
         1.3741352226e-5*N_map**4*P_map + &
         2.8567989557e-5*N_map**4*T_map**2 - &
         7.4876643038e-5*N_map**4*T_map - &
         2.0388391982e-6*N_map**4*V_map**2 + &
         4.3054466907e-5*N_map**4*V_map + &
         3.6504788687e-6*N_map**4*kappa_map**2 + &
         8.7165631487e-7*N_map**4*kappa_map + &
         1.6542743001e-5*N_map**4*mu_map**2 + &
         4.8195946039e-6*N_map**4*mu_map + &
         3.9282682647e-6*N_map**4*sigma_map**2 + &
         1.1137326431e-5*N_map**4*sigma_map + &
         2.795758112227e-5*N_map**4 + &
         1.5947545697e-6*N_map**3*P_map**2 - &
         6.9358311166e-5*N_map**3*P_map - &
         0.00014252420422*N_map**3*T_map**2 + &
         0.00039466661884*N_map**3*T_map + &
         2.15368184e-5*N_map**3*V_map**2 - &
         0.00025279065671*N_map**3*V_map + &
         4.6142483833e-6*N_map**3*kappa_map**2 - &
         2.5055687574e-5*N_map**3*kappa_map - &
         3.0424806654e-6*N_map**3*mu_map**2 - &
         4.5156027497e-5*N_map**3*mu_map - &
         1.780917608e-6*N_map**3*sigma_map**2 - &
         2.516400813e-5*N_map**3*sigma_map - &
         0.0003567127574296*N_map**3 + &
         5.9696014699e-7*N_map**2*P_map**4 - &
         1.3472490172e-5*N_map**2*P_map**3 - &
         1.0610551852e-6*N_map**2*P_map**2*T_map**2 + &
         2.0181530448e-6*N_map**2*P_map**2*T_map + &
         2.5327194907e-7*N_map**2*P_map**2*V_map**2 - &
         1.4006527233e-6*N_map**2*P_map**2*V_map + &
         5.4851851852e-7*N_map**2*P_map**2*kappa_map**2 - &
         1.320380981e-6*N_map**2*P_map**2*kappa_map + &
         1.7644666667e-7*N_map**2*P_map**2*mu_map**2 - &
         2.7894950894e-7*N_map**2*P_map**2*mu_map + &
         1.8201189815e-7*N_map**2*P_map**2*sigma_map**2 - &
         5.0510811394e-7*N_map**2*P_map**2*sigma_map - &
         6.88818634103e-6*N_map**2*P_map**2 + &
         5.0207581099e-5*N_map**2*P_map*T_map**2 - &
         0.00013814911722*N_map**2*P_map*T_map - &
         6.2792651121e-6*N_map**2*P_map*V_map**2 + &
         7.2980075931e-5*N_map**2*P_map*V_map - &
         3.7856114614e-6*N_map**2*P_map*kappa_map**2 + &
         1.2860228333e-5*N_map**2*P_map*kappa_map - &
         1.5691902399e-6*N_map**2*P_map*mu_map**2 + &
         8.2376491667e-6*N_map**2*P_map*mu_map - &
         1.3435745045e-6*N_map**2*P_map*sigma_map**2 + &
         6.0282465278e-6*N_map**2*P_map*sigma_map + &
         0.0001877522259389*N_map**2*P_map - &
         4.0442507595e-5*N_map**2*T_map**4 + &
         5.6586533058e-5*N_map**2*T_map**3 - &
         8.9548419306e-6*N_map**2*T_map**2*V_map**2 + &
         0.00014183762216*N_map**2*T_map**2*V_map + &
         1.7477041667e-7*N_map**2*T_map**2*kappa_map**2 - &
         2.2336680774e-5*N_map**2*T_map**2*kappa_map - &
         3.9516949861e-5*N_map**2*T_map**2*mu_map**2 - &
         1.428384236e-5*N_map**2*T_map**2*mu_map - &
         8.1085041667e-6*N_map**2*T_map**2*sigma_map**2 - &
         4.4004842538e-5*N_map**2*T_map**2*sigma_map + &
         0.00038258884934483*N_map**2*T_map**2 + &
         2.9970384599e-5*N_map**2*T_map*V_map**2 - &
         0.00041049796829*N_map**2*T_map*V_map + &
         6.5092115599e-6*N_map**2*T_map*kappa_map**2 + &
         3.0809800694e-5*N_map**2*T_map*kappa_map + &
         9.9551207477e-5*N_map**2*T_map*mu_map**2 + &
         1.0952167639e-5*N_map**2*T_map*mu_map + &
         2.1329980047e-5*N_map**2*T_map*sigma_map**2 + &
         8.81912525e-5*N_map**2*T_map*sigma_map - &
         0.0008911162845737*N_map**2*T_map + &
         6.9026802931e-6*N_map**2*V_map**4 - &
         4.7531336217e-5*N_map**2*V_map**3 + &
         2.5832318241e-6*N_map**2*V_map**2*kappa_map**2 - &
         1.2472907784e-6*N_map**2*V_map**2*kappa_map + &
         1.1149875079e-5*N_map**2*V_map**2*mu_map**2 - &
         2.9708960501e-6*N_map**2*V_map**2*mu_map + &
         3.2880035648e-7*N_map**2*V_map**2*sigma_map**2 + &
         7.9685785603e-6*N_map**2*V_map**2*sigma_map - &
         8.857197689645e-5*N_map**2*V_map**2 - &
         1.3905780926e-5*N_map**2*V_map*kappa_map**2 + &
         2.6425726833e-5*N_map**2*V_map*kappa_map - &
         4.4290453362e-5*N_map**2*V_map*mu_map**2 + &
         3.4602470958e-5*N_map**2*V_map*mu_map - &
         9.497372933e-6*N_map**2*V_map*sigma_map**2 - &
         8.4509070972e-6*N_map**2*V_map*sigma_map + &
         0.0007493009795633*N_map**2*V_map + &
         2.8884698866e-6*N_map**2*kappa_map**4 + &
         1.349739092e-6*N_map**2*kappa_map**3 + &
         1.7550156389e-5*N_map**2*kappa_map**2*mu_map**2 - &
         1.9786638902e-6*N_map**2*kappa_map**2*mu_map + &
         5.529520787e-6*N_map**2*kappa_map**2*sigma_map**2 + &
         1.2209966835e-5*N_map**2*kappa_map**2*sigma_map - &
         5.448370112109e-5*N_map**2*kappa_map**2 - &
         4.359847391e-5*N_map**2*kappa_map*mu_map**2 - &
         2.2737228056e-5*N_map**2*kappa_map*mu_map - &
         9.990113266e-6*N_map**2*kappa_map*sigma_map**2 - &
         5.9185131528e-5*N_map**2*kappa_map*sigma_map + &
         9.22206763018e-6*N_map**2*kappa_map + &
         3.0424263183e-5*N_map**2*mu_map**4 - &
         1.9098455668e-5*N_map**2*mu_map**3 + &
         8.54937625e-6*N_map**2*mu_map**2*sigma_map**2 - &
         4.684071842e-6*N_map**2*mu_map**2*sigma_map - &
         0.00035110649314667*N_map**2*mu_map**2 - &
         3.6261121147e-6*N_map**2*mu_map*sigma_map**2 - &
         9.2769369028e-5*N_map**2*mu_map*sigma_map + &
         0.00011212992202954*N_map**2*mu_map + &
         5.0891009441e-6*N_map**2*sigma_map**4 + &
         3.5893477645e-6*N_map**2*sigma_map**3 - &
         7.197212424173e-5*N_map**2*sigma_map**2 - &
         0.00011060069230486*N_map**2*sigma_map + &
         0.00151313669719111*N_map**2 - &
         1.1284287469e-6*N_map*P_map**4 + &
         3.0704412322e-5*N_map*P_map**3 + &
         1.6278653855e-6*N_map*P_map**2*T_map**2 - &
         3.3672619444e-6*N_map*P_map**2*T_map - &
         6.2110532065e-8*N_map*P_map**2*V_map**2 + &
         3.2172427639e-6*N_map*P_map**2*V_map - &
         2.8443321387e-7*N_map*P_map**2*kappa_map**2 + &
         7.4341916667e-7*N_map*P_map**2*kappa_map - &
         7.2252756038e-7*N_map*P_map**2*mu_map**2 + &
         8.4614527778e-7*N_map*P_map**2*mu_map - &
         1.2720654237e-7*N_map*P_map**2*sigma_map**2 + &
         2.7419097222e-7*N_map*P_map**2*sigma_map + &
         8.764064001385e-6*N_map*P_map**2 - &
         9.5804124167e-5*N_map*P_map*T_map**2 + &
         0.00027775604478*N_map*P_map*T_map + &
         1.2225588236e-5*N_map*P_map*V_map**2 - &
         0.0001716045343*N_map*P_map*V_map + &
         1.8806313889e-6*N_map*P_map*kappa_map**2 - &
         1.2701263287e-6*N_map*P_map*kappa_map + &
         1.3923449444e-5*N_map*P_map*mu_map**2 - &
         1.4186857698e-6*N_map*P_map*mu_map + &
         3.9634190278e-6*N_map*P_map*sigma_map**2 + &
         9.947051115e-6*N_map*P_map*sigma_map - &
         0.0003223815596977*N_map*P_map + &
         7.5931315992e-5*N_map*T_map**4 - &
         0.00011373913927*N_map*T_map**3 + &
         1.5275617964e-5*N_map*T_map**2*V_map**2 - &
         0.00027033311532*N_map*T_map**2*V_map - &
         9.5487976257e-6*N_map*T_map**2*kappa_map**2 + &
         6.1326942361e-5*N_map*T_map**2*kappa_map + &
         1.1264628419e-5*N_map*T_map**2*mu_map**2 + &
         0.00013788716208*N_map*T_map**2*mu_map + &
         8.4656605352e-6*N_map*T_map**2*sigma_map**2 + &
         9.7041865833e-5*N_map*T_map**2*sigma_map - &
         0.00046604057151*N_map*T_map**2 - &
         5.2633321153e-5*N_map*T_map*V_map**2 + &
         0.00082461082486*N_map*T_map*V_map + &
         1.9930343472e-5*N_map*T_map*kappa_map**2 - &
         0.00014021885993*N_map*T_map*kappa_map - &
         3.2740759028e-5*N_map*T_map*mu_map**2 - &
         0.00033633604715*N_map*T_map*mu_map - &
         2.7467490833e-5*N_map*T_map*sigma_map**2 - &
         0.0002267701814*N_map*T_map*sigma_map + &
         0.0010623078872764*N_map*T_map - &
         1.158262868e-5*N_map*V_map**4 + &
         0.00010282581987*N_map*V_map**3 + &
         2.0236346649e-7*N_map*V_map**2*kappa_map**2 - &
         7.0861126667e-6*N_map*V_map**2*kappa_map - &
         1.8571179464e-7*N_map*V_map**2*mu_map**2 - &
         2.9025647069e-5*N_map*V_map**2*mu_map - &
         2.8250510694e-6*N_map*V_map**2*sigma_map**2 - &
         1.0873061236e-5*N_map*V_map**2*sigma_map + &
         9.5035330545615e-5*N_map*V_map**2 - &
         4.7114408333e-7*N_map*V_map*kappa_map**2 + &
         3.5583995899e-5*N_map*V_map*kappa_map + &
         4.3931099764e-5*N_map*V_map*mu_map**2 + &
         9.4893199047e-5*N_map*V_map*mu_map + &
         1.9266056153e-5*N_map*V_map*sigma_map**2 + &
         8.172457216e-5*N_map*V_map*sigma_map - &
         0.00114623544365757*N_map*V_map + &
         2.5465455757e-6*N_map*kappa_map**4 - &
         1.1844938245e-5*N_map*kappa_map**3 - &
         1.2548361851e-5*N_map*kappa_map**2*mu_map**2 - &
         6.6498102778e-6*N_map*kappa_map**2*mu_map - &
         2.6413318548e-6*N_map*kappa_map**2*sigma_map**2 - &
         1.9743217083e-5*N_map*kappa_map**2*sigma_map - &
         5.237116508322e-5*N_map*kappa_map**2 + &
         2.2628180833e-5*N_map*kappa_map*mu_map**2 + &
         6.4409460169e-5*N_map*kappa_map*mu_map + &
         2.1659551389e-6*N_map*kappa_map*sigma_map**2 + &
         8.2294682962e-5*N_map*kappa_map*sigma_map + &
         0.00031141975514413*N_map*kappa_map - &
         1.1565474947e-5*N_map*mu_map**4 - &
         2.1450508636e-5*N_map*mu_map**3 + &
         6.1585477702e-6*N_map*mu_map**2*sigma_map**2 - &
         8.815663375e-5*N_map*mu_map**2*sigma_map + &
         8.443778184842e-5*N_map*mu_map**2 - &
         3.4406999306e-5*N_map*mu_map*sigma_map**2 + &
         0.00027943018423*N_map*mu_map*sigma_map + &
         0.00073132224303402*N_map*mu_map - &
         8.4378328798e-6*N_map*sigma_map**4 - &
         2.0928942447e-5*N_map*sigma_map**3 + &
         9.361717372097e-5*N_map*sigma_map**2 + &
         0.00042563590086478*N_map*sigma_map - &
         0.00207631579223133*N_map - &
         1.9577562243e-8*P_map**6 + &
         9.8981784049e-7*P_map**5 - &
         5.8363352597e-8*P_map**4*T_map**2 + &
         2.6457614122e-8*P_map**4*T_map + &
         9.0459993866e-8*P_map**4*V_map**2 + &
         2.1439092975e-7*P_map**4*V_map + &
         1.7814328446e-7*P_map**4*kappa_map**2 - &
         4.1686622901e-7*P_map**4*kappa_map + &
         5.3644855238e-8*P_map**4*mu_map**2 - &
         2.0156224591e-7*P_map**4*mu_map + &
         5.8210558734e-8*P_map**4*sigma_map**2 - &
         1.8962978248e-7*P_map**4*sigma_map + &
         4.46780827354e-7*P_map**4 - &
         1.1972281072e-6*P_map**3*T_map**2 + &
         5.3768532472e-6*P_map**3*T_map + &
         2.2417961995e-7*P_map**3*V_map**2 - &
         6.4936735747e-6*P_map**3*V_map - &
         7.4042040112e-7*P_map**3*kappa_map**2 + &
         3.1988643743e-6*P_map**3*kappa_map - &
         1.1174867493e-7*P_map**3*mu_map**2 + &
         4.9097886778e-6*P_map**3*mu_map + &
         2.5563905537e-7*P_map**3*sigma_map**2 + &
         3.0112545186e-6*P_map**3*sigma_map - &
         2.528028422697e-5*P_map**3 - &
         6.2461608542e-8*P_map**2*T_map**4 + &
         7.2160816962e-9*P_map**2*T_map**3 - &
         3.9231712963e-8*P_map**2*T_map**2*V_map**2 + &
         1.2045330835e-7*P_map**2*T_map**2*V_map - &
         1.2065046296e-7*P_map**2*T_map**2*kappa_map**2 + &
         5.5655764074e-8*P_map**2*T_map**2*kappa_map - &
         7.1769768519e-7*P_map**2*T_map**2*mu_map**2 + &
         9.214390015e-7*P_map**2*T_map**2*mu_map - &
         1.1352175926e-7*P_map**2*T_map**2*sigma_map**2 - &
         7.2727690784e-8*P_map**2*T_map**2*sigma_map + &
         9.20245283507e-7*P_map**2*T_map**2 + &
         8.0179117696e-8*P_map**2*T_map*V_map**2 + &
         4.235625e-8*P_map**2*T_map*V_map + &
         2.258923022e-7*P_map**2*T_map*kappa_map**2 - &
         4.446875e-7*P_map**2*T_map*kappa_map + &
         1.319233337e-6*P_map**2*T_map*mu_map**2 - &
         2.0462986111e-6*P_map**2*T_map*mu_map + &
         8.8850197051e-8*P_map**2*T_map*sigma_map**2 + &
         3.2751388889e-8*P_map**2*T_map*sigma_map - &
         3.083990086676e-7*P_map**2*T_map + &
         1.3373206908e-7*P_map**2*V_map**4 + &
         9.0363593389e-8*P_map**2*V_map**3 + &
         2.4463467593e-7*P_map**2*V_map**2*kappa_map**2 - &
         3.2486785978e-7*P_map**2*V_map**2*kappa_map + &
         2.7125416667e-7*P_map**2*V_map**2*mu_map**2 - &
         7.8996752457e-8*P_map**2*V_map**2*mu_map + &
         2.3869884259e-7*P_map**2*V_map**2*sigma_map**2 - &
         3.9030882886e-8*P_map**2*V_map**2*sigma_map - &
         1.676552162853e-6*P_map**2*V_map**2 - &
         3.2961961287e-7*P_map**2*V_map*kappa_map**2 + &
         1.0572459722e-6*P_map**2*V_map*kappa_map + &
         5.8846025249e-7*P_map**2*V_map*mu_map**2 - &
         1.2420694444e-7*P_map**2*V_map*mu_map + &
         1.1084042637e-7*P_map**2*V_map*sigma_map**2 + &
         5.9708680556e-7*P_map**2*V_map*sigma_map - &
         2.731802676607e-6*P_map**2*V_map + &
         1.6946466336e-7*P_map**2*kappa_map**4 - &
         1.697455568e-7*P_map**2*kappa_map**3 + &
         4.3087824074e-7*P_map**2*kappa_map**2*mu_map**2 - &
         2.6231749106e-8*P_map**2*kappa_map**2*mu_map + &
         4.1886990741e-7*P_map**2*kappa_map**2*sigma_map**2 - &
         2.1180014438e-7*P_map**2*kappa_map**2*sigma_map - &
         2.68356980142e-6*P_map**2*kappa_map**2 - &
         1.1445592205e-6*P_map**2*kappa_map*mu_map**2 + &
         2.838125e-7*P_map**2*kappa_map*mu_map - &
         7.6035506889e-7*P_map**2*kappa_map*sigma_map**2 + &
         3.2736111111e-9*P_map**2*kappa_map*sigma_map + &
         5.198700314056e-6*P_map**2*kappa_map + &
         2.570161515e-7*P_map**2*mu_map**4 - &
         4.2080255264e-7*P_map**2*mu_map**3 + &
         1.7831666667e-7*P_map**2*mu_map**2*sigma_map**2 - &
         1.3384887703e-6*P_map**2*mu_map**2*sigma_map - &
         2.364220102728e-6*P_map**2*mu_map**2 + &
         8.8640907579e-8*P_map**2*mu_map*sigma_map**2 + &
         1.2368444444e-6*P_map**2*mu_map*sigma_map + &
         4.855707265804e-6*P_map**2*mu_map + &
         9.059626753e-8*P_map**2*sigma_map**4 - &
         1.5254956604e-7*P_map**2*sigma_map**3 - &
         1.126043894964e-6*P_map**2*sigma_map**2 + &
         2.8315063203e-6*P_map**2*sigma_map - &
         2.778055205468e-6*P_map**2 - &
         2.1400747816e-6*P_map*T_map**4 + &
         4.9258105417e-6*P_map*T_map**3 + &
         5.6327936107e-7*P_map*T_map**2*V_map**2 + &
         9.3250965278e-6*P_map*T_map**2*V_map + &
         3.2710316757e-6*P_map*T_map**2*kappa_map**2 - &
         8.2942513889e-6*P_map*T_map**2*kappa_map + &
         7.9343747988e-6*P_map*T_map**2*mu_map**2 - &
         1.9181326389e-5*P_map*T_map**2*mu_map + &
         1.5642303199e-6*P_map*T_map**2*sigma_map**2 - &
         7.5503763889e-6*P_map*T_map**2*sigma_map + &
         2.871260752173e-5*P_map*T_map**2 + &
         1.4432569444e-7*P_map*T_map*V_map**2 - &
         3.6214883811e-5*P_map*T_map*V_map - &
         7.5203763889e-6*P_map*T_map*kappa_map**2 + &
         2.3241170134e-5*P_map*T_map*kappa_map - &
         1.7993693056e-5*P_map*T_map*mu_map**2 + &
         5.5177810267e-5*P_map*T_map*mu_map - &
         1.9631597222e-6*P_map*T_map*sigma_map**2 + &
         2.2474684728e-5*P_map*T_map*sigma_map - &
         0.00010697897078404*P_map*T_map + &
         4.3918577044e-7*P_map*V_map**4 - &
         6.5398867255e-6*P_map*V_map**3 - &
         1.6642529664e-6*P_map*V_map**2*kappa_map**2 + &
         4.5294820833e-6*P_map*V_map**2*kappa_map - &
         8.3413225416e-7*P_map*V_map**2*mu_map**2 + &
         3.3331961111e-6*P_map*V_map**2*mu_map - &
         5.4879252019e-7*P_map*V_map**2*sigma_map**2 + &
         2.5988245833e-6*P_map*V_map**2*sigma_map - &
         6.06167138571e-6*P_map*V_map**2 + &
         1.76840125e-6*P_map*V_map*kappa_map**2 - &
         1.2335091147e-5*P_map*V_map*kappa_map - &
         7.172135e-6*P_map*V_map*mu_map**2 - &
         1.4280271047e-5*P_map*V_map*mu_map - &
         1.7709023611e-6*P_map*V_map*sigma_map**2 - &
         1.5558439144e-5*P_map*V_map*sigma_map + &
         0.0001341572007329*P_map*V_map - &
         1.3847755834e-6*P_map*kappa_map**4 + &
         2.9663988051e-6*P_map*kappa_map**3 - &
         6.701873906e-7*P_map*kappa_map**2*mu_map**2 + &
         1.2602319444e-6*P_map*kappa_map**2*mu_map - &
         1.6684083648e-6*P_map*kappa_map**2*sigma_map**2 + &
         4.8915402778e-6*P_map*kappa_map**2*sigma_map + &
         1.830572029266e-5*P_map*kappa_map**2 + &
         7.2758208333e-6*P_map*kappa_map*mu_map**2 - &
         4.4294673496e-6*P_map*kappa_map*mu_map + &
         5.7967319444e-6*P_map*kappa_map*sigma_map**2 - &
         7.5561654674e-6*P_map*kappa_map*sigma_map - &
         6.41037679493e-5*P_map*kappa_map + &
         1.1303131108e-7*P_map*mu_map**4 + &
         4.477378242e-6*P_map*mu_map**3 - &
         8.1756005619e-8*P_map*mu_map**2*sigma_map**2 + &
         1.8019847222e-5*P_map*mu_map**2*sigma_map + &
         3.591958513889e-6*P_map*mu_map**2 + &
         3.5776194444e-6*P_map*mu_map*sigma_map**2 - &
         2.1832827595e-5*P_map*mu_map*sigma_map - &
         0.000104836086236*P_map*mu_map + &
         4.4303746065e-7*P_map*sigma_map**4 + &
         3.5105365953e-6*P_map*sigma_map**3 - &
         5.903007579301e-6*P_map*sigma_map**2 - &
         6.46457674367e-5*P_map*sigma_map + &
         0.000248152385516871*P_map + &
         9.3486231047e-7*T_map**6 - &
         1.8168892496e-6*T_map**5 - &
         1.2017117971e-7*T_map**4*V_map**2 - &
         6.0623117829e-6*T_map**4*V_map - &
         2.2138037142e-6*T_map**4*kappa_map**2 + &
         5.3907485972e-6*T_map**4*kappa_map - &
         1.5012731379e-5*T_map**4*mu_map**2 + &
         2.9320254519e-5*T_map**4*mu_map - &
         9.9591672455e-7*T_map**4*sigma_map**2 + &
         4.5407528726e-6*T_map**4*sigma_map - &
         2.1701625406048e-5*T_map**4 - &
         8.86080478e-8*T_map**3*V_map**2 + &
         1.4786332237e-5*T_map**3*V_map + &
         3.9370511477e-6*T_map**3*kappa_map**2 - &
         1.1123904654e-5*T_map**3*kappa_map + &
         1.9629299957e-5*T_map**3*mu_map**2 - &
         4.3941049344e-5*T_map**3*mu_map + &
         1.1089594981e-6*T_map**3*sigma_map**2 - &
         9.8896573442e-6*T_map**3*sigma_map + &
         4.92291899313038e-5*T_map**3 - &
         1.8208011851e-7*T_map**2*V_map**4 - &
         3.3830247075e-6*T_map**2*V_map**3 + &
         1.5875634259e-7*T_map**2*V_map**2*kappa_map**2 - &
         8.2594310191e-7*T_map**2*V_map**2*kappa_map - &
         7.0372356481e-6*T_map**2*V_map**2*mu_map**2 + &
         1.1314601854e-5*T_map**2*V_map**2*mu_map + &
         2.6109148148e-7*T_map**2*V_map**2*sigma_map**2 - &
         1.4817938429e-6*T_map**2*V_map**2*sigma_map + &
         2.584581978913e-6*T_map**2*V_map**2 + &
         8.1820206167e-6*T_map**2*V_map*kappa_map**2 - &
         2.152989125e-5*T_map**2*V_map*kappa_map + &
         3.5087114657e-5*T_map**2*V_map*mu_map**2 - &
         7.62298625e-5*T_map**2*V_map*mu_map + &
         3.7441553065e-6*T_map**2*V_map*sigma_map**2 - &
         1.9480860556e-5*T_map**2*V_map*sigma_map + &
         8.078026266135e-5*T_map**2*V_map - &
         1.4079168102e-6*T_map**2*kappa_map**4 + &
         8.494577264e-7*T_map**2*kappa_map**3 - &
         1.9940069444e-6*T_map**2*kappa_map**2*mu_map**2 - &
         3.7981316228e-6*T_map**2*kappa_map**2*mu_map + &
         1.2927453704e-7*T_map**2*kappa_map**2*sigma_map**2 - &
         5.9931564037e-6*T_map**2*kappa_map**2*sigma_map + &
         2.918584228186e-5*T_map**2*kappa_map**2 - &
         5.4461049955e-7*T_map**2*kappa_map*mu_map**2 + &
         1.7234859722e-5*T_map**2*kappa_map*mu_map - &
         1.9437451044e-6*T_map**2*kappa_map*sigma_map**2 + &
         1.7031693056e-5*T_map**2*kappa_map*sigma_map - &
         6.0924394269614e-5*T_map**2*kappa_map - &
         1.3344139356e-5*T_map**2*mu_map**4 + &
         6.9440170146e-6*T_map**2*mu_map**3 - &
         6.7769083333e-6*T_map**2*mu_map**2*sigma_map**2 + &
         3.2599345224e-6*T_map**2*mu_map**2*sigma_map + &
         0.00022442210764479*T_map**2*mu_map**2 + &
         7.6319402846e-6*T_map**2*mu_map*sigma_map**2 + &
         7.2359722222e-6*T_map**2*mu_map*sigma_map - &
         0.0003491716663951*T_map**2*mu_map - &
         1.1631111786e-6*T_map**2*sigma_map**4 + &
         2.6724248621e-6*T_map**2*sigma_map**3 + &
         1.649502809964e-5*T_map**2*sigma_map**2 - &
         6.1536219106916e-5*T_map**2*sigma_map + &
         0.000140143705596224*T_map**2 + &
         2.7736623456e-7*T_map*V_map**4 + &
         1.5654708427e-5*T_map*V_map**3 + &
         2.1280663429e-6*T_map*V_map**2*kappa_map**2 - &
         3.8522365278e-6*T_map*V_map**2*kappa_map + &
         2.1368239286e-5*T_map*V_map**2*mu_map**2 - &
         3.6494785e-5*T_map*V_map**2*mu_map + &
         2.4667129876e-8*T_map*V_map**2*sigma_map**2 + &
         8.2984694444e-7*T_map*V_map**2*sigma_map - &
         1.16151114743199e-6*T_map*V_map**2 - &
         2.4386513472e-5*T_map*V_map*kappa_map**2 + &
         7.2626637568e-5*T_map*V_map*kappa_map - &
         8.6385334444e-5*T_map*V_map*mu_map**2 + &
         0.00022158505878*T_map*V_map*mu_map - &
         6.7549275e-6*T_map*V_map*sigma_map**2 + &
         6.7690826574e-5*T_map*V_map*sigma_map - &
         0.000319030815886*T_map*V_map + &
         5.7579921563e-6*T_map*kappa_map**4 - &
         5.9128382699e-6*T_map*kappa_map**3 + &
         5.8599504703e-6*T_map*kappa_map**2*mu_map**2 + &
         1.0919901389e-5*T_map*kappa_map**2*mu_map + &
         2.3570428983e-6*T_map*kappa_map**2*sigma_map**2 + &
         1.1143115278e-5*T_map*kappa_map**2*sigma_map - &
         9.05515382865e-5*T_map*kappa_map**2 - &
         5.1905069444e-6*T_map*kappa_map*mu_map**2 - &
         4.8460056021e-5*T_map*kappa_map*mu_map - &
         1.9954263889e-6*T_map*kappa_map*sigma_map**2 - &
         4.0364821013e-5*T_map*kappa_map*sigma_map + &
         0.0002112574003288*T_map*kappa_map + &
         3.6831785874e-5*T_map*mu_map**4 - &
         2.2679927293e-5*T_map*mu_map**3 + &
         1.3165213945e-5*T_map*mu_map**2*sigma_map**2 - &
         1.5411466667e-5*T_map*mu_map**2*sigma_map - &
         0.0005305662075903*T_map*mu_map**2 - &
         1.4521058333e-5*T_map*mu_map*sigma_map**2 - &
         2.36356423e-5*T_map*mu_map*sigma_map + &
         0.0008959089906671*T_map*mu_map + &
         3.1144363024e-6*T_map*sigma_map**4 - &
         1.1900713105e-5*T_map*sigma_map**3 - &
         3.6444491206927e-5*T_map*sigma_map**2 + &
         0.000221944115643271*T_map*sigma_map - &
         0.000541037097804642*T_map + &
         1.3872452626e-7*V_map**6 + &
         2.6280928855e-6*V_map**5 + &
         1.4116924747e-6*V_map**4*kappa_map**2 - &
         2.6532439544e-6*V_map**4*kappa_map + &
         4.0826945223e-6*V_map**4*mu_map**2 - &
         6.7305962741e-6*V_map**4*mu_map + &
         3.2377252949e-7*V_map**4*sigma_map**2 - &
         2.3853637883e-7*V_map**4*sigma_map - &
         2.12902319506e-6*V_map**4 - &
         3.4146924781e-6*V_map**3*kappa_map**2 + &
         1.1308634888e-5*V_map**3*kappa_map - &
         4.3203910556e-6*V_map**3*mu_map**2 + &
         2.0522142146e-5*V_map**3*mu_map - &
         1.1130033126e-7*V_map**3*sigma_map**2 + &
         9.6641202763e-6*V_map**3*sigma_map - &
         6.9565263098929e-5*V_map**3 + &
         9.124821437e-7*V_map**2*kappa_map**4 - &
         2.374870717e-7*V_map**2*kappa_map**3 + &
         2.0551150926e-6*V_map**2*kappa_map**2*mu_map**2 + &
         2.1737964134e-6*V_map**2*kappa_map**2*mu_map + &
         1.5400502778e-6*V_map**2*kappa_map**2*sigma_map**2 + &
         2.1762389258e-6*V_map**2*kappa_map**2*sigma_map - &
         1.908307985662e-5*V_map**2*kappa_map**2 - &
         2.2644466603e-6*V_map**2*kappa_map*mu_map**2 - &
         7.2844616667e-6*V_map**2*kappa_map*mu_map - &
         1.327008481e-6*V_map**2*kappa_map*sigma_map**2 - &
         7.3647294444e-6*V_map**2*kappa_map*sigma_map + &
         3.083805354799e-5*V_map**2*kappa_map + &
         6.376029485e-6*V_map**2*mu_map**4 - &
         3.4468156835e-6*V_map**2*mu_map**3 + &
         1.9027688889e-6*V_map**2*mu_map**2*sigma_map**2 - &
         1.9749133587e-6*V_map**2*mu_map**2*sigma_map - &
         8.992990807087e-5*V_map**2*mu_map**2 - &
         1.2826181036e-6*V_map**2*mu_map*sigma_map**2 - &
         8.6458333333e-7*V_map**2*mu_map*sigma_map + &
         0.000119137833700857*V_map**2*mu_map + &
         3.1330206897e-7*V_map**2*sigma_map**4 - &
         5.9611039059e-7*V_map**2*sigma_map**3 - &
         5.61848663091e-6*V_map**2*sigma_map**2 + &
         1.0076975521136e-5*V_map**2*sigma_map - &
         2.21975659993502e-6*V_map**2 - &
         5.5542905873e-6*V_map*kappa_map**4 + &
         8.2797491074e-6*V_map*kappa_map**3 - &
         5.5680226085e-6*V_map*kappa_map**2*mu_map**2 - &
         1.0037508333e-6*V_map*kappa_map**2*mu_map - &
         5.6383374563e-6*V_map*kappa_map**2*sigma_map**2 + &
         4.7589241667e-6*V_map*kappa_map**2*sigma_map + &
         8.232677423507e-5*V_map*kappa_map**2 + &
         2.0663045e-5*V_map*kappa_map*mu_map**2 + &
         7.3770892156e-6*V_map*kappa_map*mu_map + &
         1.2777554444e-5*V_map*kappa_map*sigma_map**2 + &
         4.4264848543e-6*V_map*kappa_map*sigma_map - &
         0.0002273892174654*V_map*kappa_map - &
         1.7755825532e-5*V_map*mu_map**4 + &
         2.2503273728e-5*V_map*mu_map**3 - &
         7.4071030901e-6*V_map*mu_map**2*sigma_map**2 + &
         5.0431555833e-5*V_map*mu_map**2*sigma_map + &
         0.00023575988653791*V_map*mu_map**2 + &
         1.6386676389e-5*V_map*mu_map*sigma_map**2 - &
         4.6608524981e-5*V_map*mu_map*sigma_map - &
         0.00059482989619126*V_map*mu_map - &
         4.6406148944e-7*V_map*sigma_map**4 + &
         1.2492365373e-5*V_map*sigma_map**3 + &
         6.64618520895e-6*V_map*sigma_map**2 - &
         0.00023138236574996*V_map*sigma_map + &
         0.000737904956621858*V_map + &
         8.9855895986e-7*kappa_map**6 + &
         2.5470308117e-9*kappa_map**5 + &
         1.2059968463e-6*kappa_map**4*mu_map**2 + &
         3.1793778702e-6*kappa_map**4*mu_map + &
         1.0401467472e-6*kappa_map**4*sigma_map**2 + &
         2.1906938947e-6*kappa_map**4*sigma_map - &
         2.433005812556e-5*kappa_map**4 - &
         4.7304299279e-7*kappa_map**3*mu_map**2 - &
         6.5083289391e-6*kappa_map**3*mu_map - &
         1.4235522045e-7*kappa_map**3*sigma_map**2 - &
         5.8027556768e-6*kappa_map**3*sigma_map + &
         1.4243394357223e-5*kappa_map**3 + &
         3.9738486622e-6*kappa_map**2*mu_map**4 - &
         2.1223005863e-6*kappa_map**2*mu_map**3 + &
         1.0183000926e-5*kappa_map**2*mu_map**2*sigma_map**2 - &
         1.2959512062e-5*kappa_map**2*mu_map**2*sigma_map - &
         3.411711272594e-5*kappa_map**2*mu_map**2 - &
         1.0317703172e-5*kappa_map**2*mu_map*sigma_map**2 + &
         1.7773383333e-5*kappa_map**2*mu_map*sigma_map - &
         3.0801900036594e-5*kappa_map**2*mu_map + &
         1.5939164141e-6*kappa_map**2*sigma_map**4 + &
         1.3527431493e-6*kappa_map**2*sigma_map**3 - &
         2.054649657405e-5*kappa_map**2*sigma_map**2 - &
         2.925616248782e-5*kappa_map**2*sigma_map + &
         0.00020667647366024*kappa_map**2 - &
         7.1129482287e-6*kappa_map*mu_map**4 + &
         2.347842395e-7*kappa_map*mu_map**3 - &
         2.3796144071e-5*kappa_map*mu_map**2*sigma_map**2 + &
         1.8951372222e-5*kappa_map*mu_map**2*sigma_map + &
         4.588667312192e-5*kappa_map*mu_map**2 + &
         2.7744291667e-5*kappa_map*mu_map*sigma_map**2 - &
         4.024641369e-5*kappa_map*mu_map*sigma_map + &
         0.0001409512704825*kappa_map*mu_map - &
         2.3923579072e-6*kappa_map*sigma_map**4 - &
         6.8064892728e-6*kappa_map*sigma_map**3 + &
         2.465996737684e-5*kappa_map*sigma_map**2 + &
         0.000123292481750089*kappa_map*sigma_map - &
         0.000440486811020821*kappa_map + &
         4.9302956951e-6*mu_map**6 + &
         2.366766686e-6*mu_map**5 + &
         8.8283463137e-6*mu_map**4*sigma_map**2 - &
         1.3070610726e-5*mu_map**4*sigma_map - &
         0.0001299209556219*mu_map**4 - &
         1.7778517447e-5*mu_map**3*sigma_map**2 + &
         2.0156346188e-5*mu_map**3*sigma_map + &
         7.22110904344e-6*mu_map**3 + &
         2.7912531806e-7*mu_map**2*sigma_map**4 - &
         1.6060547415e-6*mu_map**2*sigma_map**3 - &
         4.639438308883e-5*mu_map**2*sigma_map**2 + &
         5.69529527941e-5*mu_map**2*sigma_map + &
         0.00107865381644979*mu_map**2 + &
         1.6967924396e-6*mu_map*sigma_map**4 - &
         8.6474054636e-6*mu_map*sigma_map**3 + &
         5.0242308290501e-5*mu_map*sigma_map**2 + &
         0.00011747955359353*mu_map*sigma_map - &
         0.00160943569318205*mu_map + &
         7.81942578e-7*sigma_map**6 - &
         2.6973326406e-6*sigma_map**5 - &
         1.688778013916e-5*sigma_map**4 + &
         6.520497638123e-5*sigma_map**3 + &
         9.1407247080762e-5*sigma_map**2 - &
         0.00056464306868082*sigma_map + &
         0.00111457799998979

   ! Must be greater than 0
   smax = max( rmin, smax )

end subroutine modal_smax

!=======================================================================
function lognorm_to_norm( x, lambda, zeta ) RESULT ( y )
!=======================================================================
!
! Purpose: Map a value from a lognormal distribution to standard normal.
!
!
! Author: Daniel Rothenberg, 12/02/2013
!-----------------------------------------------------------------------

   implicit NONE

!---------------------------- Arguments --------------------------------
   real(r8), intent(in) :: x
   real(r8), intent(in) :: lambda, zeta ! mean/std dev of corr. normal

   real(r8) :: y

!-----------------------------------------------------------------------
   ! In case of user stupidity...
   if ( x .le. 0._r8 ) y = -1.0_r8
   if ( x .ge. 1e6_r8 ) y = 1.0_r8

   y = ( log( x ) - lambda ) / zeta
   
end function lognorm_to_norm

!=======================================================================
function uni_to_norm( x, a, b ) RESULT ( y )
!=======================================================================
!
! Purpose: Map a value from a uniform distribution to standard normal.
!
!
! Author: Daniel Rothenberg, 12/02/2013
!-----------------------------------------------------------------------

   implicit NONE

!---------------------------- Arguments --------------------------------
   real(r8), intent(in) :: x
   real(r8), intent(in) :: a, b ! lower and upper bounds

   real(r8) :: y

!-----------------------------------------------------------------------
   if ( x .le. a ) then
     y = a
   else if ( x .ge. b ) then
     y = b
   else
     y = sqrt( 2.0_r8 )*erfinv( 2.0_r8*( x - a )/( b - a ) - 1.0_r8 )
   end if

end function uni_to_norm

!=======================================================================
function erfinv( x ) RESULT ( y )
!=======================================================================
!
! Purpose: Calculate the inverse of the error function.
!
! Utilizes the approximation developed by Mike Giles and reported in the
! book, "GPU Gems: Volume 2", chapter "Approximating the `erfinv`
! function".
!
! Author: Daniel Rothenberg, 12/02/2013
!-----------------------------------------------------------------------

   implicit NONE

!---------------------------- Arguments --------------------------------
   real(r8), intent(in) :: x

   !output
   real(r8) :: y
!------------------------ Local work space  ----------------------------
   real(r8) :: w, p
!-----------------------------------------------------------------------

   w = -1.0*log( ( 1.0_r8 - x )*( 1.0_r8 + x ) )
   if ( w < 6.250000_r8 ) then
      w = w - 3.125000_r8
      p =  -3.6444120640178196996e-21_r8
      p =   -1.685059138182016589e-19_r8 + p*w
      p =   1.2858480715256400167e-18_r8 + p*w
      p =    1.115787767802518096e-17_r8 + p*w
      p =   -1.333171662854620906e-16_r8 + p*w
      p =   2.0972767875968561637e-17_r8 + p*w
      p =   6.6376381343583238325e-15_r8 + p*w
      p =  -4.0545662729752068639e-14_r8 + p*w
      p =  -8.1519341976054721522e-14_r8 + p*w
      p =   2.6335093153082322977e-12_r8 + p*w
      p =  -1.2975133253453532498e-11_r8 + p*w
      p =  -5.4154120542946279317e-11_r8 + p*w
      p =    1.051212273321532285e-09_r8 + p*w
      p =  -4.1126339803469836976e-09_r8 + p*w
      p =  -2.9070369957882005086e-08_r8 + p*w
      p =   4.2347877827932403518e-07_r8 + p*w
      p =  -1.3654692000834678645e-06_r8 + p*w
      p =  -1.3882523362786468719e-05_r8 + p*w
      p =    0.0001867342080340571352_r8 + p*w
      p =  -0.00074070253416626697512_r8 + p*w
      p =   -0.0060336708714301490533_r8 + p*w
      p =      0.24015818242558961693_r8 + p*w
      p =       1.6536545626831027356_r8 + p*w 
   else if ( w < 16.000000_r8 ) then
      w = sqrt(w) - 3.250000_r8                          
      p =   2.2137376921775787049e-09_r8                 
      p =   9.0756561938885390979e-08_r8 + p*w           
      p =  -2.7517406297064545428e-07_r8 + p*w           
      p =   1.8239629214389227755e-08_r8 + p*w           
      p =   1.5027403968909827627e-06_r8 + p*w           
      p =   -4.013867526981545969e-06_r8 + p*w           
      p =   2.9234449089955446044e-06_r8 + p*w           
      p =   1.2475304481671778723e-05_r8 + p*w           
      p =  -4.7318229009055733981e-05_r8 + p*w           
      p =   6.8284851459573175448e-05_r8 + p*w           
      p =   2.4031110387097893999e-05_r8 + p*w           
      p =   -0.0003550375203628474796_r8 + p*w           
      p =   0.00095328937973738049703_r8 + p*w           
      p =   -0.0016882755560235047313_r8 + p*w           
      p =    0.0024914420961078508066_r8 + p*w           
      p =   -0.0037512085075692412107_r8 + p*w           
      p =     0.005370914553590063617_r8 + p*w           
      p =       1.0052589676941592334_r8 + p*w           
      p =       3.0838856104922207635_r8 + p*w                                                          
   else
      w = sqrt(w) - 5.000000_r8                          
      p =  -2.7109920616438573243e-11_r8                 
      p =  -2.5556418169965252055e-10_r8 + p*w           
      p =   1.5076572693500548083e-09_r8 + p*w           
      p =  -3.7894654401267369937e-09_r8 + p*w           
      p =   7.6157012080783393804e-09_r8 + p*w           
      p =  -1.4960026627149240478e-08_r8 + p*w           
      p =   2.9147953450901080826e-08_r8 + p*w           
      p =  -6.7711997758452339498e-08_r8 + p*w           
      p =   2.2900482228026654717e-07_r8 + p*w           
      p =  -9.9298272942317002539e-07_r8 + p*w           
      p =   4.5260625972231537039e-06_r8 + p*w           
      p =  -1.9681778105531670567e-05_r8 + p*w           
      p =   7.5995277030017761139e-05_r8 + p*w           
      p =  -0.00021503011930044477347_r8 + p*w           
      p =  -0.00013871931833623122026_r8 + p*w           
      p =       1.0103004648645343977_r8 + p*w           
      p =       4.8499064014085844221_r8 + p*w
   end if

   y = p*x   

end function erfinv
end module MITaer_activation

program main

use MITaer_activation, only: modal_smax
use shr_kind_mod, only: r8=>shr_kind_r8

real(r8) smax

call modal_smax( 1000._r8, 0.05_r8, 2._r8, 0.54_r8, 0.5_r8, &
                 283._r8, 80000._r8, smax ) 

!write (*,'(A4, 1x, F5.3, 1x, F5.1, 1x, F8.1)') "test", 0.5_r8, 283._r8, 80000._r8 
write (*,'("test (", 2(1x, F8.1), 2(1x, E8.2E2)")")') 283._r8, 80000._r8, 0.5_r8, 0.0031

write (*,*) smax
end program main