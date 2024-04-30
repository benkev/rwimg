pro formula

nx = 100
xt = -10. + 20.*dindgen(nx)/double(nx-1)

a = abs(xt)/3.0
b = dblarr(nx)
y = dblarr(nx)

for i = 0, nx-1 do begin
    if a[i] lt 1.0 then b[i] = asin(a[i]) else b[i] = !pi/2.
    if abs(a[i]) lt 0.15d then begin
        y[i] = 11./70. 
    endif else begin
        y[i] = 6.0*((b[i] - 0.50*sin(2*b[i]))/(8*a[i]^3) - (1.50*b[i] - sin(2*b[i]) + 0.1250*sin(4*b[i]))/(4*a[i]^5) $
                 +(2.50*b[i] - 1.8750*sin(2*b[i]) + 0.3750*sin(4*b[i]) - (1.0/24.0)*sin(6*b[i]))/(8*a[i]^7))

    endelse
endfor
;p = where(a lt 1.0, cntp)
;q = where(a ge 1.0, cntq)
;if cntp ne 0 then b[p] = asin(a[p])
;if cntq ne 0 then b[q] = asin(1.0)

;y = 6.0*((b - 0.5*sin(2*b))/(8*a^3) - (1.5*b - sin(2*b) + 0.125*sin(4*b))/(4*a^5) $
;       +(2.5*b - 1.875*sin(2*b)+0.375*sin(4*b) - (1.0/24.0)*sin(6*b))/(8*a^7))


window, xs=800, ys=800, /fr
plot, xt, y, bac=255, col=0, thick=3.

stop
end
