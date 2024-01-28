clear all;clf;colormap jet;

Lx = 1;
Ly = 1;
nx = 256;
ny = nx;
dt = 1e-2;
nt = 3;
pi = 3.1415926535898;


dx = Lx/nx;
dy = Ly/ny;

# Accelerated PT
tau = 1e6;
Re  = pi;
c   = 1/sqrt(2) * 0.95;
rho = Re / Lx / c * tau / dx;
th  = Lx / Re / c * tau / dx;

#tau = (min(dx,dy)^2)/4.1
#tau = a*(min(dx,dy))/4.1

x = 0.5*dx : dx : Lx-0.5*dx;
y = 0.5*dy : dy : Ly-0.5*dy;

#================ Variables and initialization
C = 0*ones(nx,ny);
% Cell loop: set params
for i = 1:nx
  for j = 1:ny
    xi = (i+0.5)*dx;
    yj = (j+0.5)*dy;
    r = sqrt( (xi-0.5)^2 + (yj-0.5)^2 );
    if r < 0.2
      C(i,j) = 1;
    endif
  endfor
endfor
Cn = C;
qx = 0*ones(nx+1, ny);
qy = 0*ones(nx, ny+1);

tic
r0 = 0;
for it = 1:nt
  printf("\n\nTime step %d, t = %d\n", it, it*dt);
  if mod(it,1) == 0
    pcolor(x,y, transpose(C)); shading interp;
    drawnow;
  endif
  Ckm1 = C;
  Ck = C;
  for i = 1:1000000

    divq = diff(qx,1,1)/dx + diff(qy,1,2)/dy;
    R = (C-Cn)./dt + divq;
    r = max(abs(R));
    r = max(r);
    if i == 2
      r0 = r
    endif
    if((r < 1e-6 || r < 1e-5*r0) && i > 1)
      break;
    endif

    Ckm1 = Ck;
    Ck = C;
    #C = C - (divq + (C-Cn)./dt) .* tau;
    C = (-R + Ck.*(rho/tau + 2*rho*th/tau/tau) - (rho*th/tau/tau) * Ckm1) ./ (rho*th/tau/tau + rho/tau);
    qx(2:end-1,:) = -diff(C,1,1)/dx;
    qy(:,2:end-1) = -diff(C,1,2)/dy;

    #qy(:,1)   = -(C(:,1) - 1)/(dy/2);
    #qy(:,end) = -(0 - C(:,end))/(dy/2);

    %qx(1,:) = 1;%h(1,:)/(dx/2);
    %qx(end,:) = -1;%(h(end,:) - 1)/(dx/2);
    if mod(i,100) == 0
      #printf(" iter %d\n", i);
      printf(" iter %d: r = %e\n", i, r);
      #pcolor(x,y, transpose(C)); shading flat; drawnow; colorbar;
    %quiver(qcx, qcy); drawnow;
    endif
  endfor

  Cn = C;
endfor
toc


