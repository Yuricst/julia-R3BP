"""
Equations of motion in R3BP
"""

using DifferentialEquations
using LinearAlgebra
using Distributed
using Printf


# -------------------------------------------------------------------------------- #
# Equations of motion for CR3BP
"""
    rhs_cr3bp_sv!(du,u,p,t)

Right-hand side expression for state-vector in CR3BP

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_cr3bp_sv!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 + z^2 );
    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # derivatives of velocities
    du[4] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[5] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    du[6] = -((1-p[1])/r1^3)*z - (p[1]/r2^3)*z;
end


"""
    rhs_cr3bp_svstm!(du,u,p,t)

Right-hand side expression for state-vector and STM in CR3BP

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_cr3bp_svstm!(du,u,p,t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 + z^2 );
    # duatives of positions
    du[1] = u[4];
    du[2] = u[5];
    du[3] = u[6];
    # duatives of velocities
    du[4] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[5] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    du[6] = -((1-p[1])/r1^3)*z - (p[1]/r2^3)*z;

    # coefficients of A matrix
    # first ~ third row
    a00, a01, a02, a03, a04, a05 = 0, 0, 0, 1, 0, 0;
    a10, a11, a12, a13, a14, a15 = 0, 0, 0, 0, 1, 0;
    a20, a21, a22, a23, a24, a25 = 0, 0, 0, 0, 0, 1;
    # fourth ~ sixth row
    a33, a34, a35 =  0, 2, 0;
    a43, a44, a45 = -2, 0, 0;
    a53, a54, a55 =  0, 0, 0;

    # STM duative
    # first row ...
    du[7]  = a00*u[7]  + a01*u[13] + a02*u[19] + a03*u[25] + a04*u[31] + a05*u[37];
    du[8]  = a00*u[8]  + a01*u[14] + a02*u[20] + a03*u[26] + a04*u[32] + a05*u[38];
    du[9]  = a00*u[9]  + a01*u[15] + a02*u[21] + a03*u[27] + a04*u[33] + a05*u[39];
    du[10] = a00*u[10] + a01*u[16] + a02*u[22] + a03*u[28] + a04*u[34] + a05*u[40];
    du[11] = a00*u[11] + a01*u[17] + a02*u[23] + a03*u[29] + a04*u[35] + a05*u[41];
    du[12] = a00*u[12] + a01*u[18] + a02*u[24] + a03*u[30] + a04*u[36] + a05*u[42];

    # second row ...
    du[13] = a10*u[7]  + a11*u[13] + a12*u[19] + a13*u[25] + a14*u[31] + a15*u[37];
    du[14] = a10*u[8]  + a11*u[14] + a12*u[20] + a13*u[26] + a14*u[32] + a15*u[38];
    du[15] = a10*u[9]  + a11*u[15] + a12*u[21] + a13*u[27] + a14*u[33] + a15*u[39];
    du[16] = a10*u[10] + a11*u[16] + a12*u[22] + a13*u[28] + a14*u[34] + a15*u[40];
    du[17] = a10*u[11] + a11*u[17] + a12*u[23] + a13*u[29] + a14*u[35] + a15*u[41];
    du[18] = a10*u[12] + a11*u[18] + a12*u[24] + a13*u[30] + a14*u[36] + a15*u[42];

    # third row ...
    du[19] = a20*u[7]  + a21*u[13] + a22*u[19] + a23*u[25] + a24*u[31] + a25*u[37];
    du[20] = a20*u[8]  + a21*u[14] + a22*u[20] + a23*u[26] + a24*u[32] + a25*u[38];
    du[21] = a20*u[9]  + a21*u[15] + a22*u[21] + a23*u[27] + a24*u[33] + a25*u[39];
    du[22] = a20*u[10] + a21*u[16] + a22*u[22] + a23*u[28] + a24*u[34] + a25*u[40];
    du[23] = a20*u[11] + a21*u[17] + a22*u[23] + a23*u[29] + a24*u[35] + a25*u[41];
    du[24] = a20*u[12] + a21*u[18] + a22*u[24] + a23*u[30] + a24*u[36] + a25*u[42];

    # fourth row ...
    du[25] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[7]  + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[13] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[19] + a33*u[25] + a34*u[31] + a35*u[37];

    du[26] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[8]  + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[14] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[20] + a33*u[26] + a34*u[32] + a35*u[38];

    du[27] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[9]  + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[15] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[21] + a33*u[27] + a34*u[33] + a35*u[39];

    du[28] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[10]  + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[16] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[22] + a33*u[28] + a34*u[34] + a35*u[40];

    du[29] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[11] + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[17] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[23] + a33*u[29] + a34*u[35] + a35*u[41];

    du[30] = ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*((u[1]+p[1])^2 *(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)^2*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[12] + ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[18] + ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[24] + a33*u[30] + a34*u[36] + a35*u[42];

    # fifth row ...
    du[31] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[7]  + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[13] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[19] + a43*u[25] + a44*u[31] + a45*u[37];

    du[32] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[8]  + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[14] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[20] + a43*u[26] + a44*u[32] + a45*u[38];

    du[33] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[9]  + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[15] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[21] + a43*u[27] + a44*u[33] + a45*u[39];

    du[34] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[10]  + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[16] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[22] + a43*u[28] + a44*u[34] + a45*u[40];

    du[35] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[11] + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[17] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[23] + a43*u[29] + a44*u[35] + a45*u[41];

    du[36] = ( 3*u[2]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[12] + ( 1 - (1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[2]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[18] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[24] + a43*u[30] + a44*u[36] + a45*u[42];

    # sixth row ...
    du[37] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[7]  + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[13] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[19] + a53*u[25] + a54*u[31] + a55*u[37];

    du[38] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[8]  + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[14] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[20] + a53*u[26] + a54*u[32] + a55*u[38];

    du[39] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[9]  + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[15] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[21] + a53*u[27] + a54*u[33] + a55*u[39];

    du[40] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[10]  + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[16] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[22] + a53*u[28] + a54*u[34] + a55*u[40];

    du[41] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[11] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[17] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[23] + a53*u[29] + a54*u[35] + a55*u[41];

    du[42] = ( 3*u[3]*((u[1]+p[1])*(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + (u[1]+p[1]-1)*p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[12] + ( 3*u[2]*u[3]*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[18] + ( -(1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^1.5 - p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^1.5 + 3*u[3]^2*((1-p[1])/( (u[1]+p[1])^2 + u[2]^2 + u[3]^2 )^2.5 + p[1]/((u[1]-1+p[1])^2 + u[2]^2 + u[3]^2)^2.5) )*u[24] + a53*u[30] + a54*u[36] + a55*u[42];
end


# -------------------------------------------------------------------------------- #
# Equations of motion for PCR3BP
"""
    rhs_pcr3bp_sv!(du,u,p,t)

Right-hand side expression for state-vector in PCR3BP

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_pcr3bp_sv!(du,u,p,t)
    # unpack state
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 );
    # derivatives of positions
    du[1] = u[3]
    du[2] = u[4]
    # derivatives of velocities
    du[3] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[4] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
end



"""
    rhs_pcr3bp_svstm!(du,u,p,t)

Right-hand side expression for state-vector and STM in PCR3BP

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_pcr3bp_svstm!(du,u,p,t)
    # unpack state
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 );
    # derivatives of positions
    du[1] = u[3]
    du[2] = u[4]
    # derivatives of velocities
    du[3] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[4] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    # second-order derivatives of pseudo-potential Uxx, Uxy, Uyy
    Uxx = 1 - (1-p[1])/r1^3 + 3*(1-p[1])*(x+p[1])^2/r1^5 - p[1]/r2^3 + 3*p[1]*(x-1+p[1])^2/r2^5;
    Uyy = 1 - (1-p[1])/r1^3 + 3*(1-p[1])*y^2/r1^5        - p[1]/r2^3 + 3*p[1]*y^2/r2^5;
    Uxy =                     3*(1-p[1])*(x+p[1])*y/r1^5 + 3*p[1]*(x-1+p[1])*y/r2^5;
    # coefficients in A matrix
    a11, a12, a13, a14 = 0,   0,    1, 0;
    a21, a22, a23, a24 = 0,   0,    0, 1;
    a31, a32, a33, a34 = Uxx, Uxy,  0, 2;
    a41, a42, a43, a44 = Uxy, Uyy, -2, 0;
    # multiply out
    du[5] = a11*u[5] + a12*u[9]  + a13*u[13] + a14*u[17];
    du[6] = a11*u[6] + a12*u[10] + a13*u[14] + a14*u[18];
    du[7] = a11*u[7] + a12*u[11] + a13*u[15] + a14*u[19];
    du[8] = a11*u[8] + a12*u[12] + a13*u[16] + a14*u[20];

    du[9]  = a21*u[5] + a22*u[9]  + a23*u[13] + a24*u[17];
    du[10] = a21*u[6] + a22*u[10] + a23*u[14] + a24*u[18];
    du[11] = a21*u[7] + a22*u[11] + a23*u[15] + a24*u[19];
    du[12] = a21*u[8] + a22*u[12] + a23*u[16] + a24*u[20];

    du[13] = a31*u[5] + a32*u[9]  + a33*u[13] + a34*u[17];
    du[14] = a31*u[6] + a32*u[10] + a33*u[14] + a34*u[18];
    du[15] = a31*u[7] + a32*u[11] + a33*u[15] + a34*u[19];
    du[16] = a31*u[8] + a32*u[12] + a33*u[16] + a34*u[20];

    du[17] = a41*u[5] + a42*u[9]  + a43*u[13] + a44*u[17];
    du[18] = a41*u[6] + a42*u[10] + a43*u[14] + a44*u[18];
    du[19] = a41*u[7] + a42*u[11] + a43*u[15] + a44*u[19];
    du[20] = a41*u[8] + a42*u[12] + a43*u[16] + a44*u[20];
end


# -------------------------------------------------------------------------------- #
# Equations of motion for ER3BP
"""
    rhs_er3bp_sv!(du,u,p,t)

ER3BP equation of motion

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = ecc, p[3] = t0
    - `t`: time
"""
function rhs_er3bp_sv!(du,u,p,t)
   # ER3BP equation of motion
   # unpack arguments
   mu = p[1]
   e  = p[2]
   t0 = p[3]
   # decompose state
   x = u[1];
   y = u[2];
   z = u[3];
   vx = u[4];
   vy = u[5];
   vz = u[6];

   # calculate radii
   r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
   r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
   # --- STATE DERIVATIVE --- #
   # position-state derivative
   du[1] = vx;
   du[2] = vy;
   du[3] = vz;

   # CR3BP pseudo-potential
   Omega_x = x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x);
   Omega_y = y - ((1-mu)/r1^3)*y - (mu/r2^3)*y;
   Omega_z = -((1-mu)/r1^3)*z - (mu/r2^3)*z;

   # ER3BP term
   ecc_factor = e*cos(t-t0) / ( 1 + e*cos(t-t0) );

   # velocity-state derivative
   du[4] =  2*vy + Omega_x - ecc_factor*Omega_x;
   du[5] = -2*vx + Omega_y - ecc_factor*Omega_y;
   du[6] =         Omega_z - ecc_factor*(Omega_z + z);
end



"""
    rhs_er3bp_svstm!(du,u,p,t)

ER3BP equation of motion with STM

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = ecc, p[3] = t0
    - `t`: time
"""
function rhs_er3bp_svstm!(du,u,p,t)
    # ER3BP equation of motion
    # unpack arguments
    mu = p[1]
    e  = p[2]
    t0 = p[3]
    # decompose state
    x = u[1];
    y = u[2];
    z = u[3];
    vx = u[4];
    vy = u[5];
    vz = u[6];
    # decompose stm
    stm = reshape(u[7:end], (6,6))';   # note Julia is column-major
    # for row = 1:6
    #     for col = 1:6
    #         stm(row,col) = u[6+col+(row-1)*6];
    #     end
    # end
    # calculate radii
    r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
    # --- STATE DERIVATIVE --- #
    # position-state derivative
    du[1] = vx;
    du[2] = vy;
    du[3] = vz;

    # CR3BP pseudo-potential
    Omega_x = x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x);
    Omega_y = y - ((1-mu)/r1^3)*y - (mu/r2^3)*y;
    Omega_z = -((1-mu)/r1^3)*z - (mu/r2^3)*z;

    # ER3BP term
    ecc_factor = e*cos(t-t0) / ( 1 + e*cos(t-t0) );

    # velocity-state derivative
    du[4] =  2*vy + Omega_x - ecc_factor*Omega_x;
    du[5] = -2*vx + Omega_y - ecc_factor*Omega_y;
    du[6] =         Omega_z - ecc_factor*(Omega_z + z);

    # --- A-MATRIX --- #
    # initialize A-matrix
    A11 = zeros(3,3);
    A12 = I(3);
    A22 = 2*[0  1 0;
             -1 0 0;
             0  0 0];
    # construct U matrix in CR3BP
    Uxx_C = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*((x+mu)^2 *(1-mu)/r1^5 + (x+mu-1)^2*mu/r2^5);
    Uyy_C = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*y^2*((1-mu)/r1^5 + mu/r2^5);
    Uzz_C = -(1-mu)/r1^3 - mu/r2^3 + 3*z^2*((1-mu)/r1^5 + mu/r2^5);
    Uxy_C = 3*y*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5);
    Uxz_C = 3*z*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5);
    Uyz_C = 3*y*z*((1-mu)/r1^5 + mu/r2^5);
    # convert to U matrix in ER3BP
    Uxx = Uxx_C - ecc_factor * (Uxx_C);
    Uyy = Uyy_C - ecc_factor * (Uyy_C);
    Uzz = Uzz_C - ecc_factor * (Uzz_C + 2);
    Uxy = Uxy_C - ecc_factor * (Uxy_C);
    Uxz = Uxz_C - ecc_factor * (Uxz_C + z);
    Uyz = Uyz_C - ecc_factor * (Uyz_C + z);

    Uderiv = [Uxx Uxy Uxz;
              Uxy Uyy Uyz;
              Uxz Uyz Uzz];
    # update A-matrix
    A = [A11    A12;
         Uderiv A22];
    # differential relation
    stmdot = A * stm;
    # store elements of A in du(7,1) ~ onwards
    du[7:end] = reshape(stmdot', (36,))
    # for row = 1:6
    #     for col = 1:6
    #         du[6+col+(row-1)*6] = stmdot(row,col);
    #     end
    # end
end



# -------------------------------------------------------------------------------- #
# Equations of motion for BCR4BP
"""
    rhs_bcr4bp_sv!(du,u,p,t)

BCR4BP equation of motion

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
    - `t`: time
"""
function rhs_bcr4bp_sv!(du,u,p,t)
    # unpack arguments
    mu, μ_3, t0, a_s, ω_s = p
    # decompose state
    x, y, z, vx, vy, vz = u

    # calculate radii
    r1 = sqrt( (x+mu)^2 + y^2 + z^2 )
    r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 )

    # sun position
    xs = a_s*cos(ω_s*t + t0)
    ys = a_s*sin(ω_s*t + t0)
    zs = 0.0
    r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x) + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs )
    du[5] = -2*vx + y - ((1-mu)/r1^3)*y      - (mu/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys )
    du[6] =           - ((1-mu)/r1^3)*z      - (mu/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs )
end



"""
    rhs_bcr4bp_svstm!(du,u,p,t)

BCR4BP equation of motion

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
    - `t`: time
"""
function rhs_bcr4bp_svstm!(du,u,p,t)
    # unpack arguments
    mu, μ_3, t0, a_s, ω_s = p
    # decompose state
    x, y, z, vx, vy, vz = u[1:6]
    stm = reshape(u[7:end], (6,6))';   # note Julia is column-major

    # calculate radii
    r1 = sqrt( (x+mu)^2 + y^2 + z^2 )
    r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 )

    # sun position
    xs = a_s*cos(ω_s*t + t0)
    ys = a_s*sin(ω_s*t + t0)
    zs = 0.0
    r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x) + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs )
    du[5] = -2*vx + y - ((1-mu)/r1^3)*y      - (mu/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys )
    du[6] =           - ((1-mu)/r1^3)*z      - (mu/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs )

    # initialize A-matrix
    A11 = zeros(3,3);
    A12 = I(3);
    A22 = 2*[0  1 0;
             -1 0 0;
             0  0 0];

    # Construct U matrix (double-deriv of potential) in CR3BP
    Uxx_C = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*((x+mu)^2 *(1-mu)/r1^5 + (x+mu-1)^2*mu/r2^5)
    Uyy_C = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*y^2*((1-mu)/r1^5 + mu/r2^5)
    Uzz_C =   - (1-mu)/r1^3 - mu/r2^3 + 3*z^2*((1-mu)/r1^5 + mu/r2^5)
    Uxy_C =                               3*y*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5)
    Uxz_C =                               3*z*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5)
    Uyz_C =                               3*y*z*((1-mu)/r1^5 + mu/r2^5)

    # construct U matrix from BCR4BP terms
    Ups_xx = 3*(μ_3/r3^5)*(x-xs)^2 - μ_3/r3^3
    Ups_yy = 3*(μ_3/r3^5)*(y-ys)^2 - μ_3/r3^3
    Ups_zz = 3*(μ_3/r3^5)*(z-zs)^2 - μ_3/r3^3
    Ups_xy = 3*(μ_3/r3^5)*(y-ys)*(x-xs)
    Ups_xz = 3*(μ_3/r3^5)*(z-zs)*(x-xs)
    Ups_yz = 3*(μ_3/r3^5)*(z-zs)*(y-ys)
    
    # combine U
    Uxx = Uxx_C + Ups_xx
    Uyy = Uyy_C + Ups_yy
    Uzz = Uzz_C + Ups_zz
    Uxy = Uxy_C + Ups_xy
    Uxz = Uxz_C + Ups_xz
    Uyz = Uyz_C + Ups_yz
    Uderiv = [Uxx Uxy Uxz;
              Uxy Uyy Uyz;
              Uxz Uyz Uzz];

    # update A-matrix
    A = [A11    A12;
         Uderiv A22];
    # differential relation
    stmdot = A * stm;
    # store elements of A in du(7,1) ~ onwards
    du[7:end] = reshape(stmdot', (36,))
end