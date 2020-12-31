"""
Function for obtaining lagrange points
"""

#using Roots

function lagrangePoint(mu::Float64)
    """Function computes lagrange points"""

    l = 1 - mu;
    # collinear points
    fl1(x) = x^5 + 2*(mu-l)*x^4 + (l^2-4*l*mu+mu^2)*x^3 + (2*mu*l*(l-mu)+mu-l)*x^2 + (mu^2*l^2+2*(l^2+mu^2))*x + (mu^3-l^3);
    xl1 = find_zero(fl1, (0,2), Bisection());

    fl2(x) = x^5 + 2*(mu-l)*x^4 + (l^2-4*l*mu+mu^2)*x^3 + (2*mu*l*(l-mu)-mu-l)*x^2 + (mu^2*l^2+2*(l^2-mu^2))*x - (mu^3+l^3);
    xl2 = find_zero(fl2, (0,2), Bisection());

    fl3(x) = x^5 + 2*(mu-l)*x^4 + (l^2-4*l*mu+mu^2)*x^3 + (2*mu*l*(l-mu)+mu+l)*x^2 + (mu^2*l^2+2*(mu^2-l^2))*x + (mu^3+l^3);
    xl3 = find_zero(fl3, (-2,0), Bisection());

    LP = zeros(5, 3);
    LP[1,1] = xl1;
    LP[2,1] = xl2;
    LP[3,1] = xl3;
    # equilateral points
    LP[4,1] = cos(pi/3)-mu;
    LP[4,2] = sin(pi/3);
    LP[5,1] = cos(pi/3)-mu;
    LP[5,2] = -sin(pi/3);
    return LP
end
