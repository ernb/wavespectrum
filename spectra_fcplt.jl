
using Plots
using SpecialFunctions

function Sgamma(freq, omegaP, l_, n_)
    omega=2pi*freq
    return omega ^ -l_ * exp(-l_ / n_ * (omega / omegaP) ^ -n_)
end
function M0_(l,n,omp)
    1/n * gamma((l-1)/n)*(l/n*omp^n)^-((l-1)/n) #verm.factor D hoort hier ook bij, maar die is apart gehouden
end  

function spectrum(Hm0_d, He10_d ::Float64)
    
 # voor numerieke integratie He10 - frequenties (40 stappen):
  freq=Array{Float64}(.01:.09/40:.1)
 # en gewichten (Simpson)
  w=Array{Int64}(ones(41))
  for i=2:2:40;w[i]=4;end
  for i=3:2:40;w[i]=2;end
  # gewenste M0 met 10% ophogen(...) omdat meting de hoogfrequente staart mist
  M0_des=1.1*Hm0_d^2 /16

 #deelspectrum zeegang resp deining, initiele piekperiodes
  Tp=[6.,6.]
  l_=[5,5]
  n_=[4,4]
  om_p = 2pi ./ Tp

  M0=M0_.(l_,n_,om_p)
  D0=M0_des ./ M0

 #stapgrootte integratie is 2pi* .09/40 ; /3 vanwege Simpson:
  Opp=Array{Float64}(undef,2)
  Opp[1]=pi*.09/20/3 * sum(w.*Sgamma.(freq,om_p[1],l_[1],n_[1]))
  Opp[2]=pi*.09/20/3 * sum(w.*Sgamma.(freq,om_p[2],l_[2],n_[2]))

 # eerst kijken of gewenste He10 uit spectra met initiele piekperioden te maken is
 #ondergrens, alles zeegang
  He10=4*sqrt(Opp[1]*D0[1])
  println("Minimale He10 = $He10")
  while He10 >= He10_d 
    Tp[1]*=0.9
    om_p = 2pi ./ Tp
    M0=M0_.(l_,n_,om_p)
    D0=M0_des ./ M0
    Opp[1]=pi*.09/20/3 * sum(w.*Sgamma.(freq,om_p[1],l_[1],n_[1]))
    He10=4*sqrt(Opp[1]*D0[1])
  end
  println("Minimale He10 = $He10")

#println("Opp = $Opp, D0 = $D0, M0 = $M0")

#bovengrens, alles deining
  He10=4*sqrt(Opp[2]*D0[2])
  println("Maximale He10 = $He10")
  while He10 <= He10_d 
    Tp[2]*=1.1
    om_p = 2pi ./ Tp
    M0=M0_.(l_,n_,om_p)
    D0=M0_des ./ M0
    Opp[2]=pi*.09/20/3 * sum(w.*Sgamma.(freq,om_p[2],l_[2],n_[2]))
    He10=4*sqrt(Opp[2]*D0[2])
  end

  println("Maximale He10 = $He10")
  println("Piekperioden: $Tp s")

# Nu lineaire combinatie van de 2 spectra bepalen: 
# zoek passende factor voor zeegang, factor deining volgt uit constant houden Hm0.
  f_zeegang=[0.,1.]
  f_s=[0.,0.]
  He10_m=100

  while abs(He10_m/He10_d -1)>.01 # op 1 cm benaderen goed genoeg

     f_s[1]=sum(f_zeegang) / 2
     # bereken overblijvende deel voor deining
     f_s[2]= 1-f_s[1]
     # He10 midpunt
     
     S_10=f_s .* D0 .* Opp
     He10_m = 4*sqrt(sum(S_10))
     #### println("Opp = $Opp, D0 = $D0, M0 = $M0, f_s = $f_s")
 
     if He10_m > He10_d
        # midpunt wordt ondergrens (hoe hoger f_zeegang hoe lager He10)
        f_zeegang[1]=f_s[1]
     else
        # midpunt wordt bovengrens
        f_zeegang[2]=f_s[1]
     end
     #### println("He10_m = $He10_m , interval $f_zeegang")
     
  end

  x=.01:.005:.5
  y1=2*pi*f_s[1]*D0[1]*Sgamma.(x,om_p[1],l_[1],n_[1]) #dichtheid per Hz is pi* per rad!
  y2=2*pi*f_s[2]*D0[2]*Sgamma.(x,om_p[2],l_[2],n_[2])
  y=y1+y2

#### println("Opp = $Opp, D0 = $D0, M0 = $M0, f_s = $f_s")
###

###
  Hm0=4*sqrt(sum(f_s .* D0 .* M0_.(l_,n_,om_p))) #dit zou nog steeds Hm0_d moeten zijn
  #Plots.GRBackend()
  He10_r=round(He10_m; digits=3)
  Hm0_r=round(Hm0; digits=3)
  plotly()

global h=plot(x,y1,
     title="Spectrum met He10=$He10_r en Hm0=$Hm0_r", 
     label="zeegang",
     xlabel="golffrequentie [Hz]",
     ylabel="energiedichtheid m2/Hz")
  plot!(x,y2,label="deining")
  plot!(x,y,label="totaal",lw=3)
display(h)
#energiebanden berekenen met Simpson, 4 stappen per bin
  w=zeros(Float64,27,109)
  for i=0:26; w[i+1,1+i*4]=1.0 ;
    w[i+1,2+i*4]=4.0 ;
    w[i+1,3+i*4]=2.0 ;
    w[i+1,4+i*4]=4.0 ;
    w[i+1,5+i*4]=1.0 ;
  end
  xx=[0.025:.0025:.295;]
  yy=zeros(Float64,length(xx))
# energiedichtheidswaarden in de steunpunten:
  for i=1:length(xx); yy[i]=2*pi*sum(f_s .* D0 .* Sgamma.(xx[i],om_p,l_,n_)); end
  Energiehoogte_cm2= w * yy .* (.0025*10000/3)
  #=
  de frequentiebins van de waarnemingen lopen tot 295 mHz, 
  dus de hieruit bepaalde M0 is te laag. De output beslaat 
  hetzelfde deel van het spectrum, waaruit nu een Hm0 bepaald wordt.
  =#
  M0ger=sum(Energiehoogte_cm2)
  Hm0ger=4*sqrt(M0ger)
  println("gerealiseerde M0_f<0.3 = $M0ger met Hm0_f<0.3 = $Hm0ger")
  #println("Ingevoerd: Hm0 = $Hm0_d m, He10 = $He10_d m. ") 
  #println("Energiehoogte_cm2: \n$Energiehoogte_cm2")
  
  return Energiehoogte_cm2

end

sp=spectrum(1.26,0.025)
const output = IOBuffer()
using REPL
const out_terminal = REPL.Terminals.TerminalBuffer(output)
const basic_repl = REPL.BasicREPL(out_terminal)
const basic_display = REPL.REPLDisplay(basic_repl)
Base.pushdisplay(basic_display)
using Plots
plotly()

display(h)
print(sp)