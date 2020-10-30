clear all
%% Initial non-calculating data

% from example 12.3 (p.752), example 12.5 (p.789) P. Kundur "Power System Stability and Control"
omega_0 = 377;
X_d = 1.81;     X_d_ht = 0.3;       X_d_2ht = 0.23;
X_q = 1.76;     X_q_ht = 0.65;      X_q_2ht = 0.25;
T_d0_ht = 8;    T_d0_2ht = 0.03;
T_q0_ht = 1;    T_q0_2ht = 0.07;
R_a = 0.003;    X_l = 0.16;
H = 3.5;        K_D = 0;
R_E = 0;        X_E = 0.65;
A_Sat = 0.031;  B_Sat = 6.93;       psi_TI = 0.8;
P = 0.9;        Q = 0.3;            E_t = 1.0;
%% Initial calculating data

% inductances
L_d = X_d;  L_d_ht = X_d_ht;    L_d_2ht = X_d_2ht;
L_q = X_q;  L_q_ht = X_q_ht;    L_q_2ht = X_q_2ht;
L_l = X_l;

% unsaturated mutual inductances:
L_ad = L_d - L_l;
L_aq = L_q - L_l;

% field winding inductance:
L_fd = L_ad*(L_d_ht - L_l)/(L_ad - L_d_ht + L_l);

% damping windings inductance:
L_1q = L_aq*(L_q_ht - L_l)/(L_aq - L_q_ht + L_l);
L_1d = -L_ad * L_fd * (L_l - L_d_2ht) / (L_ad * L_fd + L_ad * L_l - ...
    L_ad * L_d_2ht + L_fd * L_l - L_fd * L_d_2ht);
L_2q = -L_aq * L_1q * (L_l - L_q_2ht)/(L_1q * L_aq + L_1q * L_l -...
    L_1q * L_q_2ht + L_aq * L_l - L_aq * L_q_2ht);

% rotor windings inductance:
R_fd = (L_ad + L_fd) / (T_d0_ht * omega_0);
R_1d = (L_1d + (L_ad * L_fd) / (L_ad + L_fd)) / (T_d0_2ht * omega_0);
R_1q = (L_aq + L_1q) / (T_q0_ht * omega_0);
R_2q = (L_2q + (L_aq * L_1q) / (L_aq + L_1q)) / (T_q0_2ht * omega_0);

% total saturation factors K_sd, K_sq:
I_t_tild = (P + Q * 1i)' / E_t;
E_a_tild = E_t + (R_a + X_l * 1i) * I_t_tild;
psi_at = abs(E_a_tild);
if psi_at <= psi_TI
    psi_I = 0;
else
    psi_I = A_Sat * exp(B_Sat * (psi_at-psi_TI));
end
K_sd = psi_at / (psi_at + psi_I);
K_sq = K_sd;

%% Operating point parameters

% saturated values of synchronous reactances & inductances:
X_ds = K_sd * L_ad + L_l;
X_qs = K_sq * L_aq + L_l;
L_ds = X_ds;
L_qs = X_qs;

% terminal current & power factor angle:
I_t = abs(I_t_tild');
f = (angle(I_t_tild'));

% internal angle:
delta_i = rad2deg (atan((I_t * X_qs * cos((f)) -...
    I_t * R_a * sin((f))) / (E_t + I_t * R_a * cos((f))...
    + I_t * X_qs * sin((f)))));

% rotor voltages & currents:
e_d0 = E_t * sin (deg2rad(delta_i));
e_q0 = E_t * cos (deg2rad(delta_i));
i_d0 = I_t * sin (deg2rad(delta_i) + f);
i_q0 = I_t * cos (deg2rad(delta_i) + f);

% bus voltage:
E_Bd0 = e_d0 - R_E * i_d0 + X_E * i_q0;
E_Bq0 = e_q0 - R_E * i_q0 - X_E * i_d0;
E_B = sqrt(E_Bd0^2 + E_Bq0^2);

% rotor angle:
delta_0 = rad2deg(atan(E_Bd0 / E_Bq0));

% field winding current & voltage:
i_fd0 = (e_q0 + R_a * i_q0 + L_ds * i_d0) / (K_sd * L_ad);
E_fd0 = L_ad * i_fd0;

% mutual flux linkage:
psi_ad0 = K_sd * L_ad * (-i_d0 + i_fd0);
psi_aq0 = - K_sq * L_aq * i_q0;
psi_fd0 = (K_sd * L_ad + L_fd) * i_fd0 - K_sd * L_ad * i_d0;
psi_1d0 = K_sd * L_ad * (i_fd0 - i_d0);
psi_1q0 = -K_sq * L_aq * i_q0;
psi_2q0 = -K_sq * L_aq * i_q0;

%% Incremental saturation

% incremental saturation factor:
psi_at0 = psi_at;
K_sd_incr = 1 / (1 + B_Sat * A_Sat * exp(B_Sat * (psi_at0 - psi_TI)));
K_sq_incr = K_sd_incr;

% incremental saturated values of the mutual inductances:
L_ads = K_sd_incr * L_ad;
L_aqs = K_sq_incr * L_aq;

%% Calculation K_ij

L_ads_2ht = 1 / (1 / L_ads + 1 / L_fd + 1 / L_1d);
L_aqs_2ht = 1 / (1 / L_aqs + 1 / L_1q + 1 / L_2q);

k_Xd = (X_E + L_ads_2ht + L_l) / ((R_a + R_E)^2 + (X_E + L_ads_2ht + L_l)...
    * (X_E + L_aqs_2ht + L_l));
k_Xq = (X_E + L_aqs_2ht + L_l) / ((R_a + R_E)^2 + (X_E + L_ads_2ht + L_l)...
    * (X_E + L_aqs_2ht + L_l));
k_R = (R_a + R_E) / ((R_a + R_E)^2 + (X_E + L_ads_2ht + L_l)...
    * (X_E + L_aqs_2ht + L_l));

delta_0r = deg2rad(delta_0);

K_11 = - K_D / (2 * H);

K_12 = - (-E_B*L_1d*L_1q*L_2q*L_fd*(k_R^2 - k_Xd*k_Xq)*(L_ads_2ht - ...
    L_aqs_2ht)*cos(delta_0r)^2 + (2*E_B*k_R*L_1d*L_1q*L_2q*L_fd*(k_Xq +...
    k_Xd)*(L_ads_2ht - L_aqs_2ht)*sin(delta_0r) + ((psi_fd0*((k_R^2 - ...
    k_Xd*k_Xq)*L_ads_2ht + (-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + k_Xd)*...
    L_ads_2ht*L_2q + 2*k_R*L_fd*L_aqs_2ht*(L_ads_2ht*k_Xd - L_aqs_2ht*...
    k_Xd + 1/2)*psi_2q0)*L_1q + 2*k_R*psi_1q0*L_fd*L_aqs_2ht*(L_ads_2ht*...
    k_Xd - L_aqs_2ht*k_Xd + 1/2)*L_2q)*L_1d + L_1q*L_fd*((k_R^2 - k_Xd*...
    k_Xq)*L_ads_2ht + (-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + k_Xd)*L_2q*...
    psi_1d0*L_ads_2ht)*cos(delta_0r) - 2*sin(delta_0r)*(-E_B*L_1d*L_1q...
    *L_2q*L_fd*(k_R^2 - k_Xd*k_Xq)*(L_ads_2ht - L_aqs_2ht)*sin(delta_0r)/2 ...
    + ((L_ads_2ht*psi_fd0*k_R*(L_ads_2ht*k_Xq - L_aqs_2ht*k_Xq - 1/2)*L_2q ...
    - L_fd*L_aqs_2ht*((k_R^2 - k_Xd*k_Xq)*L_ads_2ht + (-k_R^2 + k_Xd*k_Xq)...
    *L_aqs_2ht - k_Xq)*psi_2q0/2)*L_1q - psi_1q0*L_fd*L_aqs_2ht*L_2q*...
    ((k_R^2 - k_Xd*k_Xq)*L_ads_2ht + (-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht - ...
    k_Xq)/2)*L_1d + L_ads_2ht*psi_1d0*k_R*L_1q*L_2q*L_fd*(L_ads_2ht*k_Xq ...
    - L_aqs_2ht*k_Xq - 1/2))) * 1/(2 * L_fd * L_1d * L_1q * L_2q * H) * E_B;

% K_12_2 * D_delta^2
K_12_2 = (-E_B*k_R*L_1d*L_1q*L_2q*L_fd*(k_Xq + k_Xd)*(L_ads_2ht - L_aqs_2ht)...
    *cos(delta_0r)^2 + (-2*E_B*L_1d*L_1q*L_2q*L_fd*(k_R^2 - k_Xd*k_Xq)*...
    (L_ads_2ht - L_aqs_2ht)*sin(delta_0r) + ((L_ads_2ht*psi_fd0*k_R*...
    (L_ads_2ht*k_Xq - L_aqs_2ht*k_Xq - 1/2)*L_2q - L_fd*L_aqs_2ht*...
    ((L_ads_2ht - L_aqs_2ht)*k_R^2 - k_Xq*(L_ads_2ht*k_Xd - L_aqs_2ht*...
    k_Xd + 1))*psi_2q0/2)*L_1q - psi_1q0*L_fd*L_aqs_2ht*L_2q*((L_ads_2ht ...
    - L_aqs_2ht)*k_R^2 - k_Xq*(L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1))/2)*...
    L_1d + L_ads_2ht*psi_1d0*k_R*L_1q*L_2q*L_fd*(L_ads_2ht*k_Xq - L_aqs_2ht...
    *k_Xq - 1/2))*cos(delta_0r) + sin(delta_0r)*(2*E_B*k_R*L_1d*L_1q*L_2q...
    *L_fd*(k_Xq + k_Xd)*(L_ads_2ht - L_aqs_2ht)*sin(delta_0r) + ((psi_fd0...
    *((L_ads_2ht - L_aqs_2ht)*k_R^2 - k_Xd*(L_ads_2ht*k_Xq - L_aqs_2ht*...
    k_Xq - 1))*L_ads_2ht*L_2q + 2*k_R*L_fd*L_aqs_2ht*(L_ads_2ht*k_Xd - ...
    L_aqs_2ht*k_Xd + 1/2)*psi_2q0)*L_1q + 2*k_R*psi_1q0*L_fd*L_aqs_2ht*...
    (L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1/2)*L_2q)*L_1d + L_1q*L_fd*...
    ((L_ads_2ht - L_aqs_2ht)*k_R^2 - k_Xd*(L_ads_2ht*k_Xq - L_aqs_2ht*...
    k_Xq - 1))*L_2q*psi_1d0*L_ads_2ht)/2)*E_B/(2*L_fd*L_1d*L_1q*L_2q*H);

K_13 = (-L_1q*L_fd*E_B*((-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + L_ads_2ht*k_R^2 ...
    - L_ads_2ht*k_Xd*k_Xq + k_Xd)*L_1d*L_2q*sin(delta_0r)/2 - k_R*L_1q*...
    L_fd*E_B*(L_ads_2ht*k_Xq - L_aqs_2ht*k_Xq - 1/2)*L_1d*L_2q*...
    cos(delta_0r) + (-(((-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + L_ads_2ht*k_R^2 ...
    - L_ads_2ht*k_Xd*k_Xq + k_Xd - k_Xq)*(L_1q*psi_2q0 + L_2q*psi_1q0)*...
    L_aqs_2ht*L_fd)/2 + L_ads_2ht*psi_fd0*k_R*L_1q*L_2q*(L_ads_2ht*k_Xq ...
    - L_aqs_2ht*k_Xq - 1))*L_1d + L_ads_2ht*psi_1d0*k_R*L_1q*L_2q*L_fd*...
    (L_ads_2ht*k_Xq - L_aqs_2ht*k_Xq - 1))*L_ads_2ht/(L_fd^2*L_1d*L_1q*L_2q*H);

% K_13_2 * D_psi_fd^2
K_13_2 = L_ads_2ht^2*(-1 + (L_ads_2ht - L_aqs_2ht)*k_Xq)*k_R/(2*L_fd^2*H);

K_14 = (-L_1q*L_fd*E_B*((-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + L_ads_2ht*k_R^2 ...
    - L_ads_2ht*k_Xd*k_Xq + k_Xd)*L_1d*L_2q*sin(delta_0r)/2 - k_R*L_1q*...
    L_fd*E_B*(L_ads_2ht*k_Xq - L_aqs_2ht*k_Xq - 1/2)*L_1d*L_2q*cos(delta_0r)...
    + (-(((-k_R^2 + k_Xd*k_Xq)*L_aqs_2ht + L_ads_2ht*k_R^2 - L_ads_2ht*...
    k_Xd*k_Xq + k_Xd - k_Xq)*(L_1q*psi_2q0 + L_2q*psi_1q0)*L_aqs_2ht*...
    L_fd)/2 + L_ads_2ht*psi_fd0*k_R*L_1q*L_2q*(L_ads_2ht*k_Xq - L_aqs_2ht...
    *k_Xq - 1))*L_1d + L_ads_2ht*psi_1d0*k_R*L_1q*L_2q*L_fd*(L_ads_2ht*...
    k_Xq - L_aqs_2ht*k_Xq - 1))*L_ads_2ht/(L_fd*L_1d^2*L_1q*L_2q*H);

% K_14_2 * D_psi_1d^2
K_14_2 = L_ads_2ht^2*(-1 + (L_ads_2ht - L_aqs_2ht)*k_Xq)*k_R/(2*L_1d^2*H);

K_15 = -(-L_1q*L_fd*E_B*L_1d*L_2q*((k_R^2 - k_Xd*k_Xq)*L_ads_2ht - ...
    L_aqs_2ht*k_R^2 + L_aqs_2ht*k_Xd*k_Xq - k_Xq)*cos(delta_0r) + 2*k_R...
    *L_1q*L_fd*E_B*L_1d*(L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1/2)*L_2q...
    *sin(delta_0r) + (((k_R^2 - k_Xd*k_Xq)*L_ads_2ht - L_aqs_2ht*k_R^2 ...
    + L_aqs_2ht*k_Xd*k_Xq + k_Xd - k_Xq)*L_ads_2ht*(L_1d*psi_fd0 + L_fd...
    *psi_1d0)*L_2q + 2*L_aqs_2ht*psi_2q0*k_R*L_1d*L_fd*(L_ads_2ht*k_Xd - ...
    L_aqs_2ht*k_Xd + 1))*L_1q + 2*L_1d*L_2q*L_aqs_2ht*L_fd*k_R*psi_1q0*...
    (L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1))*L_aqs_2ht/(2*L_fd*L_1d*L_1q^2*L_2q*H);

% K_15_2 * D_psi_1q^2
K_15_2 = -k_R*(1 + (L_ads_2ht - L_aqs_2ht)*k_Xd)*L_aqs_2ht^2/(2*L_1q^2*H);

K_16 = -(-L_1q*L_fd*E_B*L_1d*L_2q*((k_R^2 - k_Xd*k_Xq)*L_ads_2ht - L_aqs_2ht...
    *k_R^2 + L_aqs_2ht*k_Xd*k_Xq - k_Xq)*cos(delta_0r) + 2*k_R*L_1q*L_fd...
    *E_B*L_1d*(L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1/2)*L_2q*sin(delta_0r)...
    + (((k_R^2 - k_Xd*k_Xq)*L_ads_2ht - L_aqs_2ht*k_R^2 + L_aqs_2ht*k_Xd...
    *k_Xq + k_Xd - k_Xq)*L_ads_2ht*(L_1d*psi_fd0 + L_fd*psi_1d0)*L_2q + ...
    2*L_aqs_2ht*psi_2q0*k_R*L_1d*L_fd*(L_ads_2ht*k_Xd - L_aqs_2ht*k_Xd + 1))...
    *L_1q + 2*L_1d*L_2q*L_aqs_2ht*L_fd*k_R*psi_1q0*(L_ads_2ht*k_Xd -...
    L_aqs_2ht*k_Xd + 1))*L_aqs_2ht/(2*L_fd*L_1d*L_1q*L_2q^2*H);

% K_16_2 * D_psi_2d^2
K_16_2 = -k_R*(1 + (L_ads_2ht - L_aqs_2ht)*k_Xd)*L_aqs_2ht^2/(2*L_2q^2*H);

% K_123 * (D_delta * D_psi_fd)
K_123 = -L_ads_2ht*E_B*(((L_ads_2ht - L_aqs_2ht)*k_R^2 - k_Xd*(L_ads_2ht*...
    k_Xq - L_aqs_2ht*k_Xq - 1))*cos(delta_0r) - 2*k_R*(L_ads_2ht*k_Xq - ...
    L_aqs_2ht*k_Xq - 1/2)*sin(delta_0r))/(2*L_fd*H);

% K_124 * (D_delta * D_psi_1d)
K_124 = -L_ads_2ht*E_B*(((L_ads_2ht - L_aqs_2ht)*k_R^2 - k_Xd*(L_ads_2ht*...
    k_Xq - L_aqs_2ht*k_Xq - 1))*cos(delta_0r) - 2*k_R*(L_ads_2ht*k_Xq - ...
    L_aqs_2ht*k_Xq - 1/2)*sin(delta_0r))/(2*L_1d*H);

% K_125 * (D_delta * D_psi_1q)
K_125 = -(((L_ads_2ht/2 - L_aqs_2ht/2)*k_R^2 - k_Xq*(L_ads_2ht*k_Xd -...
    L_aqs_2ht*k_Xd + 1)/2)*sin(delta_0r) + k_R*cos(delta_0r)*(L_ads_2ht*...
    k_Xd - L_aqs_2ht*k_Xd + 1/2)) * L_aqs_2ht*E_B / (L_1q*H);

% K_126 * (D_delta * D_psi_2q)
K_126 = -(((L_ads_2ht/2 - L_aqs_2ht/2)*k_R^2 - k_Xq*(L_ads_2ht*k_Xd - ...
    L_aqs_2ht*k_Xd + 1)/2)*sin(delta_0r) + k_R*cos(delta_0r)*(L_ads_2ht*...
    k_Xd - L_aqs_2ht*k_Xd + 1/2))*L_aqs_2ht*E_B / (L_2q*H);

% K_134 * (D_psi_fd * D_psi_1d)
K_134 = k_R*(-1 + (L_ads_2ht - L_aqs_2ht)*k_Xq)*L_ads_2ht^2/(L_fd*L_1d*H);

% K_135 * (D_psi_fd * D_psi_1q)
K_135 = -L_ads_2ht*L_aqs_2ht*((-L_ads_2ht*k_Xd + L_aqs_2ht*k_Xd - 1)*...
    k_Xq + L_ads_2ht*k_R^2 - L_aqs_2ht*k_R^2 + k_Xd)/(2*L_1q*L_fd*H);

% K_136 * (D_psi_fd * D_psi_2q)
K_136 = -L_ads_2ht*L_aqs_2ht*((-L_ads_2ht*k_Xd + L_aqs_2ht*k_Xd - 1)...
    *k_Xq + L_ads_2ht*k_R^2 - L_aqs_2ht*k_R^2 + k_Xd)/(2*L_2q*L_fd*H);

% K_145 * (D_psi_1d * D_psi_1q)
K_145 = -L_ads_2ht*L_aqs_2ht*((-L_ads_2ht*k_Xd + L_aqs_2ht*k_Xd - 1)...
    *k_Xq + L_ads_2ht*k_R^2 - L_aqs_2ht*k_R^2 + k_Xd)/(2*L_1q*L_1d*H);

% K_146 * (D_psi_1d * D_psi_2q)
K_146 = -L_ads_2ht*L_aqs_2ht*((-L_ads_2ht*k_Xd + L_aqs_2ht*k_Xd - 1)...
    *k_Xq + L_ads_2ht*k_R^2 - L_aqs_2ht*k_R^2 + k_Xd)/(2*L_2q*L_1d*H);

% K_156 * (D_psi_1q * D_psi_2q)
K_156 = -L_aqs_2ht^2*k_R*(1 + (L_ads_2ht - L_aqs_2ht)*k_Xd)/(L_1q*L_2q*H);


K_21 = omega_0;
K_22 = 0; K_23 = 0; K_24 = 0; K_25 = 0; K_26 = 0;


K_31 = 0;
K_32 = omega_0 * R_fd * L_ads_2ht * E_B / L_fd * (k_R * cos(delta_0r) ...
    - k_Xq * sin (delta_0r));
K_32_2 = - omega_0 * R_fd * L_ads_2ht * E_B / (2 * L_fd) *...
    (k_R * sin(delta_0r) + k_Xq * cos (delta_0r));
K_32_3 = - omega_0 * R_fd * L_ads_2ht * E_B / (6 * L_fd) *...
    (k_R * cos(delta_0r) - k_Xq * sin (delta_0r));
K_33 = - omega_0 * R_fd / L_fd * (1 + (L_ads_2ht * ...
    (k_Xq * L_ads_2ht - 1)) / L_fd);
K_34 = - omega_0 * R_fd / L_fd * (L_ads_2ht * (k_Xq * L_ads_2ht - 1)) / L_1d;
K_35 = omega_0 * R_fd * k_R * L_ads_2ht * L_aqs_2ht / (L_fd * L_1q);
K_36 = omega_0 * R_fd * k_R * L_ads_2ht * L_aqs_2ht / (L_fd * L_2q);

K_41 = 0;
K_42 = L_fd * R_1d * K_32 / (R_fd * L_1d);
K_42_2 = L_fd * R_1d * K_32_2 / (R_fd * L_1d);
K_42_3 = L_fd * R_1d * K_32_3 / (R_fd * L_1d);
K_43 = R_1d / R_fd * K_34;
K_44 = - omega_0 * R_1d / L_1d * (1 + (L_ads_2ht * ...
    (k_Xq * L_ads_2ht - 1)) / L_1d);
K_45 = L_fd * R_1d * K_35 / (R_fd * L_1d);
K_46 = omega_0 * R_1d * L_ads_2ht * k_R * L_aqs_2ht/(L_1d * L_2q);

K_51 = 0;
K_52 = -omega_0 * R_1q * L_aqs_2ht * E_B / L_1q * (k_R * sin(delta_0r) ...
    + k_Xd * cos (delta_0r));
K_52_2 = -omega_0 * R_1q * L_aqs_2ht * E_B / (2*L_1q) * (k_R * cos(delta_0r) ...
    - k_Xd * sin (delta_0r));
K_52_3 = omega_0 * R_1q * L_aqs_2ht * E_B / (6*L_1q) * (k_R * sin(delta_0r) ...
    + k_Xd * cos (delta_0r));
K_53 = - R_1q * K_35 / R_fd;
K_54 = - R_1q * K_45 / R_1d;
K_55 = -omega_0 * R_1q / L_1q * (1 + (L_aqs_2ht * ...
    (k_Xd * L_aqs_2ht - 1)) / L_1q);
K_56 = -omega_0 * R_1q / L_1q * (L_aqs_2ht * (k_Xd * L_aqs_2ht - 1)) / L_2q;

K_61 = 0;
K_62 = L_1q * R_2q * K_52 / (R_1q * L_2q);
K_62_2 = L_1q * R_2q * K_52_2 / (R_1q * L_2q);
K_62_3 = L_1q * R_2q * K_52_3 / (R_1q * L_2q);
K_63 = - R_2q * K_36 / R_fd;
K_64 = - omega_0 * R_2q / L_2q * L_ads_2ht * L_aqs_2ht * k_R / L_1d;
K_65 = R_2q / R_1q * K_56;
K_66 = - omega_0 * R_2q / L_2q * (1 + (L_aqs_2ht * ...
    (k_Xd * L_aqs_2ht - 1)) / L_2q);

%% Construction of matrices A,B & N for bilinear representation

A1 = [  K_11 K_12 K_13 K_14 K_15 K_16;
        K_21 K_22 K_23 K_24 K_25 K_26;
        K_31 K_32 K_33 K_34 K_35 K_36;
        K_41 K_42 K_43 K_44 K_45 K_46;
        K_51 K_52 K_53 K_54 K_55 K_56;
        K_61 K_62 K_63 K_64 K_65 K_66];

n = size(A1,1);
n2 = n + n * n;
n3 = n2 + n * n * n;

A2 = zeros(size(A1,1),n^2);
A2 (1,:) = [0 0 0 0 0 0 ...
            0 K_12_2 K_123 K_124 K_125 K_126...
            0 K_123 K_13_2 K_134 K_135 K_136...
            0 K_124 K_134 K_14_2 K_145 K_146...
            0 K_125 K_135 K_145 K_15_2 K_156...
            0 K_126 K_136 K_146 K_156 K_16_2];
A2 (:,8) = [K_12_2; 0; K_32_2; K_42_2; K_52_2; K_62_2;];

I = eye(n, n);
II = kron(I, I);

A21 = kron(A1, I) + kron(I, A1);

B0 = zeros(size(A1,1),2);
B0(1,1) = 1/(2*H);
B0(3,2) = omega_0 * R_fd / L_ad;
B01 = zeros(size(A1,1),1);
B01(1,1) = 1/(2*H);
B201 = kron(B01, I) + kron(I, B01);
B02 = zeros(size(A1,1),1);
B02(3,1) = omega_0 * R_fd / L_ad;
B202 = kron(B02, I) + kron(I, B02);

A = zeros(n + n^2, n + n^2);
A(1:n, 1:n) = A1;
A(1:n, n+1:n2) = A2;
A(n+1:n2, n+1:n2) = A21;

B = zeros(n + n^2, 2);
B(1:n,:) = B0;

N1 = zeros(size(A));
N1(n+1:n2, 1:n) = B201;
N2 = zeros(size(A));
N2(n+1:n2, 1:n) = B202;

%% Transient process of bilinear system:

% simulation step
dt = 1e-4;

% simulation time
t = 0:dt:50.0;

% initial conditions 
x0 = zeros(size(A1,1),1);

% set the perturbation of the second state variable (rotor angle)
x0(2) = 0.5;

t_size = size(t, 2);
y_bilin = zeros(t_size, 6);
x = [x0; kron(x0,x0)];
u = zeros (2,1);

for i = 1:t_size
    y_bilin(i, :) = x(1:6);
    dx = A * x + N1 * x * u(1)+ N2 * x * u(2) + B * u;
    x = x + dx * dt;
end

figure()
plot(t, y_bilin)
title ({['Transient process of bilinear system: ']; ['\Delta\delta= ',num2str(x0(2))]});
legend('\omega_r', '\delta', '\psi_{fd}', '\psi_{1d}', '\psi_{1q}', '\psi_{2q}');
grid on