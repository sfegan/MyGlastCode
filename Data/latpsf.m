function [e,rP] = latpsf(ea,psf,class)
  P = 0.68;
  e = 10.^(2:.01:5);
%  e = 1000;

  ct=cos(60/180*pi):.005:cos(0/180*pi);
  if class == 0
    c0 = psf.psfscale(1);
    c1 = psf.psfscale(2);
    cg = psf.psfscale(5);
  else
    c0 = psf.psfscale(3);
    c1 = psf.psfscale(4);
    cg = psf.psfscale(5);
  end
  sc = sqrt((c0*(e/100).^cg).^2+c1^2);

  for ie=1:length(e)
    ei = e(ie);
    eact = interp1(0.5*(log10(ea.e(:,1))+log10(ea.e(:,2))),ea.area,log10(ei));
    eact = interp1(ea.ct(:,2),eact,ct);

    sict = interp1(0.5*(log10(psf.e(:,1))+log10(psf.e(:,2))),psf.sigma,...
	           log10(ei));
    sict = interp1(psf.ct(:,2),sict,ct);

    gact = interp1(0.5*(log10(psf.e(:,1))+log10(psf.e(:,2))),psf.gcore,...
	           log10(ei));
    gact = interp1(psf.ct(:,2),gact,ct);

    rPct = sc(ie).*sict.*sqrt(2.*gact.*((1-P).^(-1./(gact-1))-1));
    rPmax = max(rPct);
    r2 = (0:.001:100)*(rPmax^2);

    p = zeros(1,length(r2));
    for ict=1:length(gact)
      a = eact(ict);
      s = sict(ict)*sc(ie);
      g = gact(ict);
      px = (1-1/g)/(2*pi*s^2).*(1+r2/(2*g*s^2)).^(-g);
      p = p + eact(ict)*px;
    end
    N = sum(p);
    r68 = sqrt(r2(max(find(cumsum(p)<0.68*N))))*180/pi;
%    semilogy(sqrt(r2)/pi*180,p)
%    hold on

%    rPct = sc(ie).*sict.*sqrt(2.*gact.*((1-P).^(-1./(gact-1))-1))*180/pi;
%    rP(ie) = sum(eact.*rPct)./sum(eact);
    rP(ie) = r68;
  end 
