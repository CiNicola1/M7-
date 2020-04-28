$title  Pilot Implementation of MERGE with Transportation Demand

$include m7.def
set     at              All technologies /
                 hydro          Hydroelectric power (not vintaged),
                 nuc-1          Nuclear power with once-through fuel (Gen II and III),
                 nuc-adv        Generation IV nuclear power (breeder),
*mandiamo a zero subito i prossimi due , cola gas più tempo (2040)
                 coal-f         Coal-fired electric power,
                 oil-f          Oil-fired electric power,
                 gas-f          Gas-fired electric power,
                 coal-ccs       Coal-fired electric with CCS,
                 gas-ccs        Gas-fired electric with CCS,

*       Renewable electricity:
                 wind           Wind (pre-2010 vintages include all intermittent rnw)
                 solar-lc       New (moderate cost) solar -- no storage
                 biomass        New dedicated biomass generation
*questo da pompare all'inverosimile anche per fare un modello simle a una spece di storage solare 
                 solar-hc       Unlimited solar with storage (backstop)
*portare al 5%
                 oth-r          Existing other generation (geothermal - waste - not vintaged)

*       Non-electric (non-transportation in some regions) downstream direct use:
*azzeriamo il carobne ?
                 cldu           Coal for direct use -- (industrial only)
                 gsdu           Natural gas for direct use
                 lqdu           Liquid fuels for direct use
*abilitiamo ?
*                 rndu           Non-liquid Renewable fuels for direct use
                 neb-hc         Non-electric backstop

*       Passenger transportation vehicles
                 evhcl          Electric vehicles (for regions without transport module)

                 icev           internal combustion engine vehicle
                 phev           plug-in hybrid electric vehicle
                 elcv           full electric vehicle
                 cngv           compressed natural gas vehicle
*questo lo possiamo aumentare insieme all elettirico 
                 bksv           backstop (e.g. H2) vehicle/

        et(at)          Electricity technologies
                        /hydro, oth-r,
                         wind, solar-lc, biomass, solar-hc,
                         nuc-1, nuc-adv,
                         coal-f, oil-f, gas-f,
                         coal-ccs, gas-ccs/

        nt(at)          Nonelectric energy supply technologies
                        /cldu,lqdu,gsdu,evhcl,neb-hc/

        pvt(at)         Passenger Vehicle Technologies
                        /icev,phev,elcv,cngv,bksv/

*       Upstream fossil fuel production:

        rsc              Resources
                        /coal  Coal supply (no resource depletion accounting)
                         oil   Oil supply (conventional only),
                         gas   Natural gas supply
                         ur    Uranium
                         /

        fos(rsc)         Fossil fuel resources
                         /coal, oil, gas/

        x(rsc)           Exhaustible resources
                         /oil, gas, ur/
*cosa facciamo ? ipotizziamo risorse infinite in europa ?
*tipo gas nel mediterraneo , carbone germania polonia , petrolio mare del nord
*        xf(x)            Exhaustible fossil fuel resources
*                         /oil, gas/

*        Transportation fuels

        tf              Transportation fuels /
                         elec    Electricity
                         lqtr    Liquid transportation fuels
                         cng     Compressed natural gas
                         bstr    Backstop transportation fuel (e.g. H2)/

        ntf(tf)         Non-electric transport fuels /lqtr, cng, bstr/,

        lf              Liquid fuels
                        /petro  Petroleum-based liquids,
                         bfuel  Biomass-based liquids,
                         synf   Coal-based synthetic liquids/,

        nlf(lf)         Non-petroleum liquid fuels /bfuel,synf/,

        ghg             Greenhouse gases /co2,ch4,n2o,slf,llf/,

        ogg(ghg)        Other Greenhouse gases /ch4,n2o,slf,llf/,

        fgas(ghg)       /slf,llf/,

        ogs             Non-energy OGG sectors (read from nco2data.gdx),

        time            All time periods (starting year)
                        /1975,1980,1985,1990,1995,2000,2005,
                              2010,2020,2030,2040,2050,2060,2070,2080,2090,
                         2100,2110,2120,2130,2140,2150,2160,2170,2180,2190,2200/,

        htp(time)       Historic time periods
                        /1975,1980,1985,1990,1995,2000,2005/,
*dobbiamo dire a che anno partiamo , e se si riesce con uno step di 5 anni
        ctp(time)       Climate model time periods
                        /
                         1990,
                         2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,
                         2100,2110,2120,2130,2140,2150,2160,2170,2180,2190,2200/,

        cpp(ctp)        Climate model projection periods

        ecp(ctp)        Projection periods for both economic and climate model

        ctfix(ctp)      Fixed years in climate model
                        /
                         1990
$if not set ninety       2000
                        /

        tp(ctp)         Economic model time periods
                        /
                         1990
                         2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100/,

        hp(tp)          Economic model time periods with observed data
                        /1990, 2000/,

        pp(tp)          Economic model projection periods

        tfix(tp)        Time periods for which economic variables are fixed
                        /
                         1990
$if not set ninety       2000
$if %tenfx%==yes         2010
$if set twentyfx         2020
$if set thirtyfx         2020
$if set thirtyfx         2030
                         /,

        tbase(tp)       Base year for economic model (dynamically assigned),

        tlast(tp)       Last economic time period /2100/;

parameter       yr(ctp)        First year of time period
                nyper(ctp)     Number of years per period;

yr(ctp) = 1990 + 10 * (ord(ctp) - 1);
nyper(ctp) = 10;

* Economic model projection periods are all non-fixed periods:
pp(tp) = yes$(not tfix(tp));
* Base year for economic model is latest fixed period:
tbase(tp) = yes$(yr(tp) eq smax(tfix,yr(tfix)));


alias(ctp,cttp);

set     p0(tp)          First endogenous year;

*       This assigns p0 to be 2020 when 2010 is fixed:

        loop(tp, p0(tp)=yes$(pp(tp) and (yr(tp) le smin(pp, yr(pp)))););


set     iter            Negishi iterations  /it0*%iter%/,

        xlc(et)                Coal-fired electricity technologies subject to expansion limits
                        /coal-f, coal-ccs/

        xlg(et)                Gas-fired electricity technologies subject to expansion limits
                        /gas-f, gas-ccs/

        dle(et)                Electricity technologies subject to decline limits
                        /coal-f, coal-ccs, gas-f, gas-ccs, wind, biomass, solar-lc, solar-hc/

        nuc(et)                Nuclear electricity technologies
                        /nuc-1, nuc-adv/

        ccs(et)                Electricity technologies with CCS
                        /coal-ccs, gas-ccs/

        rnw(et)         Non-hydro non-waste renewable electricity technologies
                        /wind, solar-lc, biomass, solar-hc/

        irnw(et)        Intermittent non-hydro renewables
                        /wind, solar-lc/
*togliere o alzare questi limiti 
        crnw(et)        Constrained non-hydro renewables (all but backstop and oth)
                        /wind, solar-lc, biomass/
*metterli nel modello , o almeno togliere da qui solar hc , solar e nuc avanzato 
        etn(et)         New Electricity technologies (no deployment in 2010)
                        /nuc-adv, coal-ccs, gas-ccs, solar-lc, biomass, solar-hc/

        etf(et)         Fossil-fuel Electricity technologies
                        /coal-f, oil-f, gas-f, coal-ccs, gas-ccs/

        xersc(et,rsc)   Mapping of resources to electric technologies
                        /coal-f.coal, coal-ccs.coal,
                         gas-f.gas, gas-ccs.gas,
                         oil-f.oil, nuc-1.ur, nuc-adv.ur/
*mettere il coal 
        xln(nt)         Nonelectric technologies subject to expansion limits
                        /neb-hc/

        dln(nt)         Nonelectric technologies subject to decline limits
                        /cldu, neb-hc/

        xnfos(nt,fos)   Mapping of fossil resources to non-electric technologies
                        /cldu.coal, gsdu.gas/

        xlfos(lf,fos)   Fossil-based liquid fuels
                        /petro.oil, synf.coal/

        xlt(pvt)        Passenger vehicle technologies subject to expansion limits
                        /phev,elcv,cngv,bksv/

        trd             Tradeable goods
                        /nmr, crude, ngpl, lng, crt-ceq, crt-co2, crt-ch4, crt-n2o, crt-slf, crt-llf/

        xtrd(x,trd)     Trade of exhaustible resources
                        /oil.crude, gas.ngpl, gas.lng/

        xghg(ghg,trd)   Trade of emissions rights
                        /co2.crt-co2, ch4.crt-ch4, n2o.crt-n2o, slf.crt-slf, llf.crt-llf/

        abx             OGG abatement cost levels
                        /ab0*ab4/

        sec             Energy sectors
                        /elec,nele/

* Vintaging structure for electric sector

        v               Technology vintages
                        /1990,2000,2010,2020,2030,2040,"2050+"/,

        xlv(v)          Vintages that cannot expand (all but 2050+)
                        /1990,2000,2010,2020,2030,2040/,
        dlv(v)          Vintages subject to decline limits (all but 1990)
                        /2000,2010,2020,2030,2040,"2050+"/,

        vnw(v)          New technology vintages
        vex(v)          Existing technology vintages

        etv(et)         Electricity Technologies with vintaging
                        /wind, solar-lc, biomass, solar-hc,
                         nuc-1, nuc-adv,
                         coal-f, oil-f, gas-f,
                         coal-ccs, gas-ccs/,

        vtp(v,tp)       Time periods during which vintage v is active

        vnuc(v,tp)    Time periods in which nuclear vintage is active
;

parameter        vyr(v)          First year of vintage
                 nuclife         Lifetime of nuclear plant /60/;

vyr(v)  = sum(sameas(v,tp), yr(tp));
vyr("1990") = 1990;
vyr("2050+") = 2050;

loop(tbase, vnw(v)$(vyr(v) > yr(tbase)) = yes;
            vex(v)$(vyr(v) <= yr(tbase)) = yes;
);

vtp(v,tp)$(yr(tp) >= vyr(v)) = yes;
vnuc(v,tp)$((yr(tp) - vyr(v)) < nuclife) = yes;

* Multiple states of world for endogenous uncertainty analysis

set        sw              State of world /ref/
           st(* ,* )       Combinations of sw and tp requiring constraints,
           swm(*,*,*)      Set membership;

alias (sw,ssw,sow);
*facciamo scomparire i francesi seeeeehhhhhh
PARAMETER       prob(*)        Probability of each state of world /ref 1/,
                sp(*,*)        Previous state of world;

st(pp,sw) = yes;
swm(sw,sw,ctp) = yes;
sp(sw,tp) = 0;

alias (tp,ttp);
alias (lf,llf);

* Read the data tables:

$include readdata

$onempty
set
$if not set out                    out(tp,rg) time periods-regions outside coalition / /
$if not exist pfix_%pfix%.gms      pfix(tp,rg) time periods-regions fixed to reference case / /
$if not exist psub_%psub%.gms      psub(tp,rg) time periods-regions with full subsidy / /
$if not exist taxout_%taxout%.gms  taxout(tp,rg) time periods-regions with zero ghg tax / /
$if not exist qtyout_%qtyout%.gms  qtyout(tp,rg) time periods-regions with no qty limit and no trade / /
                                   tnbc(tp,rg) time periods over which banking-borrowing may occur / /
;
$offempty

$if set out $include out_%out%.gms

$if exist pfix_%pfix%.gms $include pfix_%pfix%.gms
$if %pfix%==yes pfix(tp,rg)$out(tp,rg) = yes;

$if exist psub_%psub%.gms $include psub_%psub%.gms
$if %psub%==yes psub(tp,rg)$out(tp,rg) = yes;

$if exist taxout_%taxout%.gms $include taxout_%taxout%.gms
$if %taxout%==yes taxout(tp,rg)$out(tp,rg) = yes;

$if exist qtyout_%qtyout%.gms $include qtyout_%qtyout%.gms
$if %qtyout%==yes qtyout(tp,rg)$out(tp,rg) = yes;


* * * * * * * * * * * * * * * * CALIBRATION * * * * * * * * * * * * * * * * *

parameter
        y0(tp,rg)       Benchmark Output
        i0(tp,rg)       Benchmark investment
        c0(tp,rg)       Benchmark consumption
        k0(tp,rg)       Benchmark capital
        e0(tp,rg)       Benchmark electric energy
        n0(tp,rg)       Benchmark non-electric energy
        t0(tp,rg)       Benchmark passenger vehicle service
        tc0(tp,rg)      Benchmark cost of transportation service
        nt0(tp,rg)      Benchmark passenger vehicle non-electric energy use

        rho(tp,rg)      Primal elasticity parameter
        kgdp(rg)        Initial CAPITAL-GDP RATIO
        depr(rg)        Annual percent depreciation
        spda(rg)        Speed of adjustment
        speed(tp,rg)    Period adjustment speed
        pnref(tp,rg)            Reference price of non-electric energy - dollars per GJ
        peref(tp,rg)            Reference price of electric energy - mills per KWH
        esub(tp,rg)             Elasticity between K-L and E-N
        kpvs(rg)                Capital value share parameter
        elvs(tp,rg)             Electric value share parameter
        theta(tp,rg)            Energy value share parameter

        xntax(rg)               Existing tax rate on nonelectric energy
        xetax(rg)               Existing tax rate on electric energy

        decf(tp,rg)             Max decline factor (annual) for PE-DL -  SYNF -  NE-BAK
        nshf(rg)                Max market share factor for EX and NX technologies
        expf(rg)                Max expansion factor (annual) for XLE technologies
        nxpf(rg)                Max expansion factor (annual) for XLN technologies

        ogpd(rg)                Oil-gas domestic delivery price differential - dollars per GJ
        clgdp(rg)               Coal-GDP growth elasticity

        dfactcurr(ctp,rg)       Current annual utility discount factor
        udf(ctp,rg)             Utility discount factor for period TP
        L(tp,rg)                Current labor force (efficiency units)
        ln(tp,rg)               Labor force - new vintage
        kl(tp,rg)               Capital-labor index

        yref(tp,rg)             Reference value of gross output
        kref(tp,rg)             Reference value of capital
        eref(tp,rg)             Reference quantity of electricity - TKWH
        nref(tp,rg)             Reference quantity of nonelectric energy

        ynref(tp,rg)            Reference path of new vintage output
        knref(tp,rg)            Reference path of new vintage capital
        enref(tp,rg)            Reference path of new vintage electric
        nnref(tp,rg)            Reference path of new vintage non-electric

        aeeifac(sec,tp,rg)       AEEI reduction factor

        pvpi(trd,tp,sw)         Present value prices of tradeables
        pvpt(tp,rg,sw)          Present value prices of transportation
        nwt(rg)                 Negishi weights
        nwtitr(iter,rg)         Negishi weights at each iteration
        taxrev(tp,rg,sw)        Energy tax revenues for each region
        taxitr(tp,rg,iter,sw)   Tax revenues at each iteration
        iprev(tp,rg,sw)         Insecurity price premium revenue recycling
        ghgtax(ghg,tp,rg,sw)    GHG tax - dollars per ton of gas
        ghgitr(ghg,tp,rg,iter,sw)  GHG tax at each iteration
        ghgsbd(ghg,tp,rg,sw)    Subsidy payment for non-participants (always > 0)
        nashsbd(ghg,tp,rg,sw)   Subsidy payment for Nash equilibrium (always > 0)
        kterm(rg,sw)            Terminal capital stock
;
kgdp(rg)  = macro("kgdp", rg);
depr(rg)  = macro("depr", rg)/100;
pnref(tp,rg) = macro("pnref", rg);
esub(tp,rg) = macro("esub", rg);

kpvs(rg)  = macro("kpvs", rg);
elvs(tp,rg)  = macro("elvs", rg);
xntax(rg) = macro("xntax", rg);
xetax(rg) = macro("xetax", rg);
decf(tp,rg)  = macro("decf", rg);
decf(tp,rg)$(yr(tp)>=2040) = 0.95;
nshf(rg)  = macro("nshf",rg);
expf(rg)  = macro("expf", rg);
nxpf(rg)  = macro("nxpf", rg);
clgdp(rg) = macro("clgdp", rg);
ogpd(rg)  = macro("ogpd", rg);

* Initialize Negishi weights
loop(tbase, nwt(rg)   = gdp(tbase,rg) / sum(rrg, gdp(tbase,rrg)););
display nwt;

k0(tfix,rg)  =  kgdp(rg)* gdp(tfix,rg);
spda(rg)      =  1 - depr(rg);
speed(tp,rg)  =  spda(rg)**10;

rho(tp,rg) = (esub(tp,rg) - 1)/esub(tp,rg);
taxrev(tp,rg,sw) = 0;
iprev(tp,rg,sw) = 0;
ghgsbd(ghg,tp,rg,sw) = 0;
nashsbd(ghg,tp,rg,sw) = 0;

aeeifac(sec,tbase,rg) = 1;
loop(pp(tp), aeeifac(sec,tp,rg) = ((1 - aeei(tp-1,rg))**10) * aeeifac(sec,tp-1,rg));
display aeeifac;

* * * * * * * * * * * * * Assigning fixed year values * * * * * * * * * * *

loop(tfix(tp),

     e0(tp,rg) = sum(et, pe0(et,tp,rg));

     n0(tp,rg) = sum(nt, pn0(nt,tp,rg));

     t0(tp,trans(rg)) = sum(pvt, pt0(pvt,tp,rg));

*      Assume no electricity used by passenger vehicles in base year(s)

     nt0(tp,trans(rg)) = sum((ntf,pvt),pt0(pvt,tp,rg) * epvkt_us(tp,ntf,pvt) * epvkt_r(pvt,rg));

*      Calculate initial investment to satisfy steady-state condition

     i0(tp,rg) = gdp(tp,rg) * kgdp(rg) * (grow(tp,rg) + depr(rg));

* Costs of transportation are subtracted from consumption
* Energy costs are assumed to represent 10% of total costs

     tc0(tp,trans(rg)) = 1.11*sum(pvt,pt0(pvt,tp,rg)*vcst(pvt,tp,rg))/1000;

     c0(tp,rg) = gdp(tp,rg) - i0(tp,rg) - tc0(tp,rg);

     y0(tp,rg) = i0(tp,rg) + c0(tp,rg)
                    + ( sum(et, pe0(et,tp,rg)*ecst0(et,rg))
                      + sum(nt, pn0(nt,tp,rg)*ncst(nt,tp,rg))
                      + sum(x, ex0(x,tp,rg)*xcst0(x,rg))
                      + sum(x, ra0(x,tp,rg)*acst0(x,rg))
                      - sum(trd, ntx0(trd,tp,rg) *  macro("intpr",rg))
* Add only non-energy costs for transport, since energy costs for base year
* are included in oil and gas market
                      + sum(pvt,pt0(pvt,tp,rg)*vcst(pvt,tp,rg))$trans(rg)
                    )/1000;
);

*       Adjust elvs to account for transport-induced changes in the
*       share of nonelectric energy in the base year:

loop(tbase,
         elvs(tp,trans(rg)) = ((n0(tbase,rg) + nt0(tbase,rg))/n0(tbase,rg)) * (elvs(tp,rg)/(1-elvs(tp,rg))) /
                           (1 + ((n0(tbase,rg) + nt0(tbase,rg))/n0(tbase,rg)) * (elvs(tp,rg)/(1-elvs(tp,rg))));

);

* Establish benchmark quantities for production function calibration

loop(tbase,
* Labor input (fixed) follows benchmark GDP path, indexed to base year
      L(tp,rg)$(tbase(tp) or pp(tp))    = gdp(tp,rg) / gdp(tbase,rg);
* Benchmark capital grows at same rate
      kref(tp,rg)$(tbase(tp) or pp(tp)) = L(tp,rg) * k0(tbase,rg);
* Benchmark energy grows at same rate less AEEI adjustments
      eref(tp,rg)$(tbase(tp) or pp(tp)) = L(tp,rg) * e0(tbase,rg) * aeeifac("elec",tp,rg);
      nref(tp,rg)$(tbase(tp) or pp(tp)) = L(tp,rg) * n0(tbase,rg) * aeeifac("nele",tp,rg);
);

* In simplified two-region world with transportation module turned off,
* we can specify a target energy path and sequentially recalibrate to it:
$if set recal parameter e_trg(tp,rg), n_trg(tp,rg);
$if set recal execute_load '%outdir%%recal%_trg.gdx', e_trg, n_trg;
$if set recal eref(tp,rg)$(tbase(tp) or pp(tp)) = e_trg(tp,rg);
$if set recal nref(tp,rg)$(tbase(tp) or pp(tp)) = n_trg(tp,rg);

* Calculate reference price of electricity (stationary if nref/eref is stationary):
peref(tp,rg)$(tbase(tp) or pp(tp)) =  pnref(tp,rg) * (elvs(tp,rg)/(1-elvs(tp,rg))) * (nref(tp,rg)/eref(tp,rg));

* Reference output assumes reference energy prices
yref(tp,rg)$(tbase(tp) or pp(tp)) =  gdp(tp,rg) +
                                    (peref(tp,rg)*eref(tp,rg) +  pnref(tp,rg)*nref(tp,rg))/1000;
* Define new vintage labor input
loop(pp(tp),
        LN(tp,rg)    = L(tp,rg) - L(tp-1,rg)*speed(tp-1,rg);

* Production function can be calibrated to gross or new vintage benchmarks
* Using gross values is isomorphic to ETA-MACRO
* To use new vintage values, derive from gross values:
        knref(tp,rg) = kref(tp,rg) - kref(tp-1,rg)*speed(tp-1,rg);
        enref(tp,rg) = eref(tp,rg) - eref(tp-1,rg)*speed(tp-1,rg);
        nnref(tp,rg) = nref(tp,rg) - nref(tp-1,rg)*speed(tp-1,rg);
        ynref(tp,rg) = yref(tp,rg) - yref(tp-1,rg)*speed(tp-1,rg);
);

* Energy value share parameter depends on gross vs. new vintage calibration:

$if not set gross $set gross yes
$if set newvint $set gross no
$if not set newvint $set newvint no


theta(pp,rg) =
         ((peref(pp,rg)*eref(pp,rg) +  pnref(pp,rg)*nref(pp,rg))/(1000*yref(pp,rg)))$%gross% +
         ((peref(pp,rg)*enref(pp,rg) +  pnref(pp,rg)*nnref(pp,rg))/(1000*ynref(pp,rg)))$%newvint%;

* When applying a policy to a sequentially recalibrated baseline (as in 1990), we must
* read the benchmark values and share parameters from the baseline run:

$if set rclbase execute_load '%outdir%%rclbase%_out.gdx', yref, ynref, kref, knref, eref, enref, nref, nnref, peref, pnref, theta, elvs;

*       Define value shares and elasticity parameters for transportation
*       demand module:

parameter       uref(tp,rg)     Reference final demand,
                alpha(tp,rg)    Value share of non-transport,
                cref(tp,rg)     Reference consumption,
                tref(tp,rg)     Reference transport demand (trillion vkt),
                tcref(tp,rg)    Reference transportation costs (trillion $),
                esubt           Elasticity between C and T /%esubt%/,
                deprt           Putty-clay depreciation rate /0.1/,
                rhot            Transport demand elasticity,
                speedt          Adjustment speed of transport demand;


* Benchmark transportation costs
tcref(tp,rg) = 0;
tref(tp,trans(rg)) = tdem(tp,rg);
loop(tbase, tcref(tp,trans(rg)) = (tdem(tp,rg) / t0(tbase,rg)) * tc0(tbase,rg););

* Old way to calculate reference consumption:  cref(tp,rg) = c0(rg) * L(tp,rg);
* With non-stationary growth rate, investment share changes over time, so
* total consumption does not grow at same rate as gdp.
* In transport model, cref refers to consumption net of transport costs.

cref(tp,rg) = gdp(tp,rg) - kref(tp,rg)*(grow(tp,rg) + depr(rg)) - tcref(tp,rg) ;

alpha(tp,tpclay(rg)) = cref(tp,rg)/(tcref(tp,rg) + cref(tp,rg));

display alpha;

uref(tp,tpclay(rg)) = cref(tp,rg)/alpha(tp,rg);
speedt  =  (1-deprt)**10;
rhot = (esubt - 1)/esubt;

* Show breakdown of competing claims equation in the reference path:
parameter ccref(tp,rg,*);
ccref(tp,rg,"con") = cref(tp,rg);
ccref(tp,rg,"inv") = gdp(tp,rg) - cref(tp,rg);
ccref(tp,rg,"ec") = peref(tp,rg)*eref(tp,rg) +  pnref(tp,rg)*nref(tp,rg);

* Intertemporal optimality condition is applied to calculate each period's
* annual utility discount factor individually.

grow(ctp,rg)$(not tp(ctp)) = 0.01;
mpc(ctp)$(not tp(ctp)) = mpc("2100");
dfactcurr(ctp,rg)  =  (1+grow(ctp,rg))/(1+mpc(ctp));

udf(tbase,rg)  =  1;
loop(tp$(not tfix(tp)),
      udf(tp,rg)  =  udf(tp-1,rg) * dfactcurr(tp-1,rg)**10;
);
loop(ctp$(not tp(ctp)),
      udf(ctp,rg)  =  udf(ctp-1,rg) * dfactcurr(ctp-1,rg)**10;
);


* Scale discount factors so that they sum to unity for each region
* over 2200 horizon (scaled back below to 2100 if damages are turned off):

udf(ctp,rg)  =  udf(ctp,rg)/(sum(tp$(not tfix(tp)), udf(tp,rg)) + sum(cttp$(not tp(cttp)),udf(cttp,rg)));
udf(tfix,rg) = 0;

* * * * * * * * * * * Prepare technology data * * * * * * * * * *
*---------------------------------------------------------------------------------------------------------------------
* Electric Sector Parameters
*qui sifa sul serio
*--------------------------------------------------------------------------------------------------------------
parameter
        ecst(et,tp,rg)        Levelized electricity cost for non-vintaged technologies
*modificare le tabelle
        ecst_v(et,v,tp,rg)    Levelized electricity cost for vintaged technologies
        cece_v(et,v,tp,rg)    Carbon emissions coefficient for vintaged technologies;

* Use observations for new 2000 vintage heat rate (and new 2010 vintage in case of oil)
elec_htrt(etv,tech,"2000",rg) = htrt_v(etv,"2000","2000",rg);
elec_htrt("oil-f",tech,"2010",rg) = htrt_v("oil-f","2010","2010",rg);

* For all existing vintages, fix heat rate in projection periods to base year level
loop(tbase, htrt_v(etv,vex,pp,rg) = htrt_v(etv,vex,tbase,rg););

* For all existing vintages, cost excludes capital component
loop(tbase, ecst_v(etv,vex,pp,rg) = ecst0(etv,rg););

* Costs for non-vintaged technology
ecst(et,tp,rg)$(not etv(et)) = ecst0(et,rg);

*       Specify new technology scenario
*si toglie qui il carbone ?
set     scenario(et,tech)       Describes technology scenario;
scenario("coal-f","%fossil%") = yes;
scenario("gas-f","%fossil%") = yes;
scenario("oil-f","%fossil%") = yes;
*should we let expanding not iv gen nuc ?
scenario("nuc-1","%nuc%") = yes;
scenario("nuc-adv","%nuc%") = yes;
scenario(ccs,"%ccs%") = yes;
scenario(rnw,"%rnw%") = yes;

* Costs for new technologies are given in 2006 dollars; convert to 2000 dollars
parameter deflator  Adjustment from 2006 dollars to 2000 dollars /0.8572/;

* Make sure heat rates are monotonically decreasing - fix later
* loop(vtp(vnw(v+1),tp+1), htrt_v(etf,v+1,pp,rg) = min(htrt_v(etf,v,pp,rg) +

* Create new vintages based on performance characteristics of new capacity

loop(tbase,
htrt_v(etv,vtp(vnw,pp),rg) = sum(tech$scenario(etv,tech), elec_htrt(etv,tech,vnw,rg));
ecst_v(etv,vtp(vnw,pp),rg) = deflator * (
* Non-energy cost
                                         sum(tech$scenario(etv,tech), ecst_non(etv,tech,vnw,rg)) +
* Fuel cost (zero if included elsewhere)
                                         htrt_v(etv,vnw,pp,rg) * (ecst_fuel(etv,pp,rg) +
* Storage cost (zero for non CCS technologies)
                                         scst(etv,pp,rg) * cprt(etv,vnw,rg) * sum(xersc(etv,fos), cec(fos))) );
);

* Calculate carbon emissions coefficients for vintaged technologies
cece_v(etv,vtp(v,pp),rg) = htrt_v(etv,v,pp,rg) * (1 - cprt(etv,v,rg)) * sum(xersc(etv,fos), cec(fos));
* To avoid double counting, set cece("oil-f") to 0, since all oil is tracked through PR
cece_v("oil-f",v,tp,rg) = 0;
*-----------------------------------------------------------------------------------------------------------
*qui si fa la svolta 
* Set heat rates for non-fossil technologies to 10:
htrt_v(et,v,tp,rg)$(not etf(et)) = 10;

* Upper bounds for new additions of individual technologies
*here we must change 
parameter        ecap(et,tp,rg)         Upper bounds on new vintage capacity;

ecap(etv,pp,rg)   = sum(tech$scenario(etv,tech), elec_cap(etv,tech,pp,rg));

* Upper bounds for groups of technologies
*lasciamo?(è usa)
parameter       rnwusacap(tp)           ETAC Renewable Cap for USA
                nucusacap(tp)           ETAC Nuclear Cap for USA;

* Implement ETAC's renewable limits for the USA
*   2020 limit of 50 GWe installed capacity for solar and wind combined
*   With 0.35 capacity factor, this equals 50 * 8760 * 0.35 ~ 0.15 TKWh
*   2030 limit is 70 GWe ~ 0.21 TKWh
* For earlier time steps, assume that actual upper bound was 50% higher than observed levels

rnwusacap(tp) = inf;
rnwusacap(tp)$(yr(tp) < 2020) = 1.5 * sum(usa(rg), pe0("wind",tp,rg));
rnwusacap("2020") = 0.15;
rnwusacap("2030") = 0.21;

* In limited pessimistic case, total nuclear in US is limited to 2000 production level

nucusacap(tp) = inf;
$if set nuclim nucusacap(tp) = pe0("nuc-1","2000","usa");

set    rsk      Set of items subject to societal risk  /n,c,p/;
parameter
        rsktax(rsk, tp,rg,sw)        Implicit tax based on risk premium for nuclear CCS and plutonium
        wtp0(rsk, rg)                Referece non-market cost of nuclear and CCS,
        refshr(rsk, rg)              Reference share of nuclear and CCS;

*       EMF22 runs used wtp0 = 10.  ETAC runs initially used same value,
*       but lowered it to 5 to compensate for increased market costs as a sensitivity.

*       Calibrate external cost of nuclear power to base year
*       nuclear share in the Annex-B countries.  Assume that
*       wtp in non-annexb is lower initially but will be idential
*       to that of the US when nuclear share reaches equivalent
*       levels.  For CCS, all regions are symmetric at 10%.

wtp0("n", anb) = 10;
wtp0("n", fsu) = 5;
wtp0("n", nanb) = 5;
wtp0("c",rg) = wtp0("n",rg);
wtp0("p",rg) = wtp0("n",rg);

refshr("n",rg) = sum(tbase, pe0("nuc-1",tbase,rg)/sum(et,pe0(et,tbase,rg)));
refshr("n",nanb) = pe_hs("nuc-1","2005","usa") / sum(et,pe_hs(et,"2005","usa"));
*refshr("n",rg)$(refshr("n",rg) eq 0) = pe_hs("nuc-1","2000","usa") / sum(et,pe_hs(et,"2000","usa"));
refshr("c",rg) = 0.2;
refshr("p",rg) = 0.05;

* Only cldu has non-electric capacity constraint in the future
*limiti sul carbone per riscladamento 
parameter        ncap(nt,tp,rg)         Upper bounds on non-electric technologies;

ncap(nt,tbase,rg) = pn0(nt,tbase,rg);
ncap(nt,pp,rg)    = inf;

loop(pp(tp+1),    ncap("cldu",tp+1,rg) =
            ncap("cldu",tp,rg)*(1+grow(tp,rg)*clgdp(rg))**10);
ncap("cldu",tp,rg)$(yr(tp)>=2050)  = ncap("cldu","2050",rg);

rcap(lf,tbase,rg) = pr0(lf,tbase,rg);

DISPLAY ncap;
*---------------------------------------------------------------------------------------------------------------------------
* Transportation Sector Parameters

parameter        epvkt(tf,pvt,tp,rg)    Energy requirement per VKT (GJ per 1000 km)
                 tfcst(tf,tp,rg)        Transportation fuel costs ($ per GJ);

* Improvement in vehicle efficiency is specified in epvkt_us (independent of AEEI):
epvkt(tf,pvt,tp,rg) = epvkt_us(tp,tf,pvt) * epvkt_r(pvt,rg);

tfcst("bstr",tp,rg) = ncst("neb-hc",tp,rg);

* Allow increase in base vehicle cost as incomes rise:

vcst("icev",tp,rg)$(ppc(tp,rg) < 20) = (vcst_z("icev","2000","usa") + sum(xzn(zn,rg), vcst_z("icev","2000",zn)))/2 -
                                       (vcst_z("icev","2000","usa") - sum(xzn(zn,rg), vcst_z("icev","2000",zn)))/2 *
                                        cos((ppc(tp,rg)/20)*pi);
vcst("icev",tp,rg)$(ppc(tp,rg) ge 20) = vcst_z("icev",tp,"usa");

*       Define electric vehicle assumptions for regions in which we
*       do not include passenger vehicles explicitly:

parameter       evlimit(tp)     Upper bound on EV share of nonelectric energy
        / 2010  0.0003, 2020 0.009, 2030 0.03, 2040 0.06, 2050 0.09, 2060 0.135 /;
loop(tp$((evlimit(tp)>0) and (evlimit(tp+1)=0)), evlimit(tp+1) = evlimit(tp););

$if set noev evlimit(tp)=0;

parameter       eveff(tp,rg)    Electric vehicle efficiency (for regions without PV submodel);

*       Rough calculation of the electric energy requirement (Tkwh) per
*       EJ of non-electric energy replaced:

*       eveff = 0.3e-12 Tkwh/mile * 30 mile/gallon * 1/1.3e-10 gallon/EJ
eveff(pp,rg) = 9/130;
*---------------------------------------------------------------------------------------------------------------------------
* Other technology/macroeconomic parameters


parameter cstexp(trd)  Unit export cost / nmr 0, crude 0.33, ngpl 0.5, lng 2.0/
          ipprem(trd,tp,rg)      Insecurity price premium
          lngx(rg)       Expansion factor for LNG imports and exports
          lngf(rg)       Fixed fraction of gas market for first year of LNG;

loop(ghg, loop(xghg(ghg,trd), cstexp(trd) = 0.1;););
ipprem("crude",tp,rg) = 2;

*  These parameters govern expansion of LNG when pipeline option is turned on

lngx(rg) = expf(rg);
lngf(rg) = 0.1;

*  If pipeline option is turned off, LNG expansion constraints are removed
*  All imports are in LNG category and follow standard gasx/gasm bounds

*  If pipeline option is turned on, gasx/gasm bounds are removed
$if set ngpl gasx(tp,rg) = +inf; gasm(tp,rg) = +inf;
*-------------------------------------------------------------------------------------------------------------------
* Define non-energy / non-CO2 emissions

parameter cempcp(tp,rg) Cement demand per capita (tons);

cempcp("1990",rg) = cement("1990",rg) / pop("1990",rg);
cempcp("2000",rg) = cement("2000",rg) / pop("2000",rg);
cempcp("2010",rg) = cement("2010",rg) / pop("2010",rg);
* Assume long-run average of 0.5 tons per person, higher for low-income regions
loop(tp$(yr(tp) >= 2020), cempcp(tp,rg) = 0.5 * cempcp(tp-1,rg) + 0.5 *
                                 ( (0.5)$(ppc(tp,rg) > 20) +
                                   (1.0)$(ppc(tp,rg) le 20) );
);
cement(tp,rg) = pop(tp,rg) * cempcp(tp,rg);

* Capture technology may be used to produce cement at a high cost:
parameter cemccs_cost(tp,rg) Cost per ton of carbon captured from cement;

cemccs_cost(tp,rg) = 1000;

* Demand for international bunker fuels grows, but emissions can be abated
parameter bunker(tp)     International bunker fuel emissions -- billion tons C
          bunker_min(tp) Minimum bunker fuel CO2 (after abatement) -- billion tons C
          bunker_abcst(tp) Cost per ton C of bunker fuel abatement;

* From ORNL data, growth rate is approximately 0.6 the global GDP growth rate since 2000)

bunker("1990") = 0.132;
bunker("2000") = 0.219;
loop(tp$(yr(tp) ge 2010), bunker(tp+1) = bunker(tp) * (1 + 0.6 * ((sum(rg, gdp(tp+1,rg))/sum(rg, gdp(tp,rg)))**0.1 - 1))**nyper(tp););

* High cost ability to abate bunker fuel emissions:

bunker_min(tp) = bunker(tp);
loop(pp(tp), bunker_min(tp+2) = 0.9 * bunker_min(tp+1););
bunker_abcst(tp) = 3666;

* Baseline emissions for OGGs are extrapolated linearly from 2020 through 2100; no
* change thereafter. Assumes equal time intervals (decades).

loop(ctp,
     oggbline(ghg,ogs,tp,rg)$(yr(tp) ge 2030) =
        oggbline(ghg,ogs,tp-1,rg) + 0.5*(oggbline(ghg,ogs,"2020",rg) - oggbline(ghg,ogs,"2000",rg));
     oggbline(ghg,ogs,tp,rg) =  max(0, oggbline(ghg,ogs,tp,rg))
    );

* N2O Baseline stays constant
oggbline("n2o",ogs,tp,rg)$(yr(tp) ge 2030) = oggbline("n2o",ogs,"2020",rg);

* Abatement is scaled with baseline emissions and adjusted for technical change
oggabt(ogg,ogs,abx,tp,rg)$((yr(tp) ge 2030) and oggbline(ogg,ogs,"2020",rg)) =
         abmlt(tp) * oggabt(ogg,ogs,abx,"2020",rg) * oggbline(ogg,ogs,tp,rg)/oggbline(ogg,ogs,"2020",rg);
oggabt(ogg,ogs,abx,tp,rg) = min(oggabt(ogg,ogs,abx,tp,rg), oggbline(ogg,ogs,tp,rg));

* Define marginal abatement limits for each cost class
* Zero cost reductions are included in the first cost category
ablim(ghg,"ab1",tp,rg) =  sum(ogs, oggabt(ghg,ogs,"ab1",tp,rg));
ablim(ghg,abx+2,tp,rg) =  sum(ogs, oggabt(ghg,ogs,abx+2,tp,rg) - oggabt(ghg,ogs,abx+1,tp,rg));
ablim(ghg,abx,tp,rg) =  max(0, ablim(ghg,abx,tp,rg));

* Define minimum possible emissions (baseline - maximal abatement) in OGG (native gas units):
oggmin(ogg,ogs,tp,rg) = (oggbline(ogg,ogs,tp,rg) - smax(abx, oggabt(ogg,ogs,abx,tp,rg)))/(1000*gwp(ogg));

* Compute actual coefficients for energy-related CH4 emissions in 1990 and 2000 (proportional to consumption).
* Assume world drops to half of 2000 U.S. level by 2050 (EPA's abatement curves not used):

ch4rate(fos,tp,rg)$(yr(tp) le 2000 and foscon0(fos,tp,rg))  =  ch4en(fos,tp,rg)/foscon0(fos,tp,rg);
ch4rate(fos,tp,rg)$(yr(tp) ge 2050) = 0.5 * sum(rrg, ch4rate(fos,"2000",rrg)$(ord(rrg) eq 1));
loop(tp$(yr(tp) < 2050 and yr(tp) > 2000),
         ch4rate(fos,tp,rg)  =  ch4rate(fos,tp-1,rg) - 0.2 * (ch4rate(fos,"2000",rg) -
                                                              ch4rate(fos,"2050",rg)) );
mlev0(tfix,rg) = sum(fos, foscon0(fos,tfix,rg)*ch4rate(fos,tfix,rg))/1000;

parameter sh(ctp,rg)            shares of global carbon emissions rights;
$ontext
The 2000 shares are proportional to 2000 regional emissions.
The 2050 shares (and thereafter) are proportional to 2000 populations.
The shares between 2000 and 2050 are determined by linear interpolation.
Shares sum to unity in each period.
$offtext

sh("2000",rg)   =     clev0("2000",rg)/sum(rrg,clev0("2000",rrg));
sh("2050",rg)   =     pop("2000",rg)/sum(rrg,pop("2000",rrg));

sh("2010",rg)   =  .8*sh("2000",rg)  +  .2*sh("2050",rg);
sh("2020",rg)   =  .6*sh("2000",rg)  +  .4*sh("2050",rg);
sh("2030",rg)   =  .4*sh("2000",rg)  +  .6*sh("2050",rg);
sh("2040",rg)   =  .2*sh("2000",rg)  +  .8*sh("2050",rg);

sh(ctp,rg)$(yr(ctp) > 2050) = sh("2050",rg);
display sh;
*-----------------------------------------------------------------------------------------------------------
* Initialize xshr for one-world resources

parameter xshr(rsc,tp,rg,sw)     Regional share of exhaustible resource consumption;

loop(tbase,
  xshr(fos,tp,rg,sw) = foscon0(fos,tbase,rg) / sum(rrg, foscon0(fos,tbase,rrg));
  xshr(fos,tp,rg,sw) = pe0("nuc-1",tbase,rg) / sum(rrg, pe0("nuc-1",tbase,rrg));
);

*  Initialize ghgtax: used to model a price-based policy, or to
*  to collect subsidy adjustments for second-best policies

ghgtax(ghg,tp,rg,sw) = 0;

* Load price policy from a previous run if specified:

$if set tax parameter totemit_tax(ghg,tp,sw), trdbal_tax(trd,tp,sw);
$if set tax execute_load '%outdir%%tax%_out.gdx', totemit_tax=totemit.m, trdbal_tax=trdbal.m;
$if set tax ghgtax(ghg,pp,rg,sw)$(not taxout(pp,rg)) = -1000 * totemit_tax(ghg,pp,sw)/abs(trdbal_tax("nmr",pp,sw));

display ghgtax;
$if set ceq ghgtax(ogg,tp,rg,sw) = gwp(ogg) * ghgtax("co2",tp,rg,sw);

* This parameter can be used to model a quantity-based policy
parameter carlim(tp,rg)                Annual carbon limits  - billion tons C
          ogglim(ghg,tp,rg)     Annual OGG limits  - billion tons of gas
          ceqlim(tp,rg)                Annual limits in billion tons of carbon-eq (using GWPs);

carlim(tp,rg) = +inf;
ogglim(ghg,tp,rg) = +inf;

* Load quantity policy from a previous run if specified:

$if set qty parameter clev_qty(tp,rg,sw), ccem_qty(tp,rg,sw), othcabt_qty(tp,sw), mlev_qty(tp,rg,sw), othemit_qty(ghg,tp,rg,sw);
$if set qty execute_load '%outdir%%qty%_out.gdx', clev_qty=CLEV.L, ccem_qty=CCEM.L, othcabt_qty=OTHCABT.L mlev_qty=MLEV.L, othemit_qty=OTHEMIT.L;
$if set qty carlim(pp(tp),rg)$(not qtyout(tp,rg)) = clev_qty(tp,rg,"ref") + ccem_qty(tp,rg,"ref");
$if set qty ogglim(ghg,pp(tp),rg)$(not qtyout(tp,rg)) = mlev_qty(tp,rg,"ref")$sameas(ghg,"ch4") + othemit_qty(ghg,tp,rg,"ref");

$if set frankel $include frankel.gms

ceqlim(tp,rg) = carlim(tp,rg) + sum(ogg(ghg), ogglim(ghg,tp,rg) * gwp(ghg));

$if set ceq carlim(tp,rg) = +inf; ogglim(ghg,tp,rg) = +inf;
$if not set ceq ceqlim(tp,rg) = +inf;

*       If we are running the EMF22 scenarios for the US, read
*       the assumptions here:

$if set emf22us $include emf22us.gms

* Adjust units for ogglim constraint to increase tolerance
parameter        ghgtol(ghg);
ghgtol(ogg) = 1;
ghgtol("llf") = 1000;

* If mu, tau, or zeta is set, model will calculate minimum feasible threshold:

$if set mu parameter mu RF penalty /%mu%/ ;
$if not set mu parameter mu /0/;
$if set tau parameter tau ATP penalty /%tau%/ ;
$if not set tau parameter tau /0/;
$if set zeta parameter zeta CON penalty /%zeta%/ ;
$if not set zeta parameter zeta /0/;


* * * * * * * * * * * * * * * VEQ * * * * * * * * * * * * * * *

*  This section of the code contains the variable and equation definitions.
*  Upper case letters denote decision variables.  Except for stock variables,
*  these represent annual rates.  Unit: trillion 2000 dollars unless indicated
*  otherwise.

POSITIVE VARIABLES
        Y(rg,tp,sw)     Production
        C(rg,ctp,sw)    Consumption - trillion dollars
        I(rg,tp,sw)     Investment - trillion dollars
        K(rg,tp,sw)     Capital stock
        E(rg,tp,sw)     Electric energy
        N(rg,tp,sw)     Non-electric energy
        T(rg,tp,sw)     Passenger transportation services

        U(rg,ctp,sw)    Utility (for model with passenger transport)

        YN(rg,tp,sw)    New production
        CN(rg,tp,sw)    New vintage goods consumption
        KN(rg,tp,sw)    New capital stock
        EN(rg,tp,sw)    New electric energy
        NN(rg,tp,sw)    New non-electric energy
        UN(rg,tp,sw)    New vintage utility
        TN(rg,tp,sw)    New-vintage transportation

        PE(et,*,rg,sw)          Production of electric energy (TKWh)
        PEV(et,v,*,rg,sw)       Vintages of electric energy (TKWh)
        PN(nt,tp,rg,sw)         Production (use) of nonelectric energy services (EJ)
        PT(pvt,tp,rg,sw)        Production of pass. trans. services - vehicle km travelled
        PR(lf,tp,rg,sw)         Production of refined liquids (EJ)

        EX(x,tp,rgw,sw)         Extraction (production) of fossil fuels (EJ)
        CX(x,tp,rgw,sw)         Cumulative production as fraction of total resource base
        XCST(x,tp,rgw,sw)       Extraction cost per unit
        RA(x,tp,rgw,sw)         Reserve additions
        ACST(x,tp,rgw,sw)       Reserve addition cost per unit
        CA(x,tp,rgw,sw)         Cumulative reserve additions as a fraction of undiscovered resources
        URSC(x,tp,rgw,sw)       Undiscovered resources
        PRSV(x,tp,rgw,sw)       Proven reserves

        CG(tp,sw)               Cumulative global consumption of coal
        SC(tp,rg,sw)            Cumulative geologically sequestered carbon from CCS technologies (billion tons C)

        CLEV(tp,rg,sw)          Energy-related Carbon emissions - billion tons C
        CCEM(tp,rg,sw)          Cement-related Carbon emissions - billion tons C
        CEMCCS(tp,rg,sw)        Cement production using capture technology -- billion tons cement
        MLEV(tp,rg,sw)          Energy-related methane emissions - billion tons CH4
        CEQ(tp,rg,sw)           Total emissions in billion tons C-eq
        EM(ghg,ctp,sw)          World energy-related emissions (CO2 and CH4 only) - billion tons C or CH4
        OTHEMIT(ghg,tp,rg,sw)   Non-energy-related GHG emissions - billion tons of gas
        ABATE(ghg,abx,tp,rg,sw) OGG abatement - million tce
        OTHCABT(tp,sw)          Non-energy CO2 abatement (billion tons C)
        TOTEM(ghg,ctp,sw)       total global anthropogenic emissions (energy and non-energy less abatement)
*** where we can go carbon negative 
*        AFF(ctp,ps,sw)          afforestation scenario - fraction adopted
*        CRLX(ctp,rg,sw)         carbon limit relaxation - billion tons
*        ORLX(ogg,ctp,rg,sw)     OGG limit relaxation - million tce
*-------------------------------------------------------------------------------------------------------------------------------
*-------------------------------------------------------------------------------------------------------------------------------
*here is the GOD 
        EXPRT(trd,tp,rg,sw)     Exports
        IMPRT(trd,tp,rg,sw)     Imports

        TPE(rsc,tp,rg,sw)       Total primary energy consumption (EJ)

        RFMAX                   Radiative Forcing limit (for diagnostic testing)
        TPMAX                   Temperature limit (for diagnostic testing)
        CONMAX                  Concentration limit (for diagnostic testing)
        ELF(rg,ctp,sw)          Economic loss factor
;


VARIABLES
        NWEL                    Negishi welfare
        EC(rg,tp,sw)            Energy cost - trillion dollars
        TC(rg,tp,sw)            Passenger transport costs - trillion dollars
        AC(rg,ctp,sw)           Abatement costs (set to zero for reference case)
        MD(rg,ctp,sw)           Market damages - trillion dollars
        NONTAX(rg,tp,sw)        Non-energy related tax payments
        NETLAND(tp,sw)          Net land use change emissions
*        AFFEMIT(tp,sw)          (Negative) Global emissions associated with afforestation - billion tons C
        NTX(trd,tp,rg,sw)       Net exports
        NBC(tp,rg,sw)           Net banking of carbon emissions rights
        CBC(tp,rg,sw)           Cumulative banked credits;

EQUATIONS

*       MACRO submodel:

        nweldf                  Negishi welfare definition
        newcap(rg,tp,sw)        New capital
        newprod(rg,tp,sw)       New production
        newelec(rg,tp,sw)       New electric energy
        newnon(rg,tp,sw)        New non-electric energy
        totalcap(rg,tp,sw)      Total capital stock
        totalprod(rg,tp,sw)     Total production
        termk(rg,tp,sw)         Terminal condition on investment and capital stock

*       Putty-clay transport demand:

        newdemand(rg,tp,sw)     New vintage final demand
        totalutl(rg,tp,sw)      Total utility
        newcon(rg,tp,sw)        New non-energy consumption
        newtrn(rg,tp,sw)        New transport

*       ETA submodel:

        suptrans(rg,tp,sw)      Supply of transportation services
        supliq(rg,tp,sw)        Supply of refined liquid fuels (for transport)
        supvint(et,rg,tp,sw)    Vintage supply
        supelec(*,*,sw)         Supply of electricity
        supnon(*,*,sw)          Supply of non-electric energy
        concoal(rg,tp,sw)       Consumption of coal
        conoil(rg,tp,sw)        Consumption of oil
        congas(rg,tp,sw)        Consumption of gas
        connuc(rg,tp,sw)        Consumption of uranium
        supx(x,rg,tp,sw)        Supply of oil and gas
        wsupx(x,tp,sw)          World supply of oil and gas
        cxdef(x,rgw,tp,sw)      Cumulative production
        xcstdef(x,rgw,tp,sw)    Extraction cost (rises with cumulative production)
        cadef(x,rgw,tp,sw)      Cumulative discoveries
        acstdef(x,rgw,tp,sw)    Discovery cost (rises with cumulative production)
        cgdef(tp,sw)            Cumulative global conusmption of coal
        scdef(rg,tp,sw)         Cumulative stored carbon
        egfrac(rg,tp,sw)        Gas fraction of electric energy
        ecfrac(rg,tp,sw)        Coal fraction of electric energy
        ccsfrac(rg,tp,sw)       Ccs - fraction of electric energy supplied by CCS tech
        intfrac(rg,tp,sw)       Intermittent share of electricity
        constrnw(rg,tp,sw)      Constrained renewable share of electricity
        rnwusa(tp,sw)           Upper bound on renewables in 2020 and 2030 in US
        evfrac(rg,tp,sw)        Electric vehicle share of non-electric energy
        nucusa(tp,sw)           Upper bound on nuclear in US in pessimistic case
        nucusa20(sw)            Limit on 2010-20 nuclear builds in US (ETAC)
        gfrac(rg,tp,sw)         Gas fraction of nonelectric energy
*        lffrac(lf,rg,tp,sw)     Limit on non-petroleum liquid fuel shares
        costnrg(rg,tp,sw)       Cost of energy
        costtrn(rg,tp,sw)       Cost of passenger transportation services
        costabt(rg,tp,sw)       Cost of non-energy abatement activities
        nontaxdef(rg,tp,sw)     Sum of payments for non-energy GHG emissions taxes
        cc(rg,tp,sw)            Capacity constraint
        trdbal(trd,tp,sw)       Global trade balance
        xptdef(trd,rg,tp,sw)    Triggers definition of positive exports
        iptdef(trd,rg,tp,sw)    Triggers definition of positive imports
        pplxpt(rg,tp,sw)        Pipeline exports
        pplipt(rg,tp,sw)        Pipeline imports
        dece(rg,dle,tp,sw)      Decline rate of dle technologies
        decev(rg,et,v,tp,sw)    Decline rate of technology vintages
        fosdec(et,v,rg,tp,sw)   Decline of existing fossil generation (forced)
        nucdec(et,v,rg,tp,sw)   New nuclear monotonicity constraints
        decn(rg,dln,tp,sw)      Decline rate of dln technologies
        dect(rg,pvt,tp,sw)      Decline rate of passenger vehicle technologies
        decr(rg,lf,tp,sw)       Decline rate of refined liquid fuels
        decx(rg,x,tp,sw)        Decline rate of fossil fuel extraction
        expev(rg,et,v,tp,sw)    Expanstion rate of technology vintages
        expnucanb(et,rg,tp,sw)  Expansion of Nuclear in Annex B
        expnucnanb(et,rg,tp,sw) Expansion of Nuclear in Non Annex B
        expc(rg,tp,sw)          Expansion rate of xlc technologies
        expg(rg,tp,sw)          Expansion rate of xlg technologies
        exprnw(et,rg,tp,sw)     Expansion rate of rnw technologies
        expccs(rg,tp,sw)        Expansion rate of CCS technologies
        expn(rg,xln,tp,sw)      Expansion rate of xln technologies
        expt(rg,pvt,tp,sw)      Expansion rate of pvt technologies
        expr(rg,lf,tp,sw)       Expansion rate of refined liquid fuels
        expcm(rg,tp,sw)         Expansion of capture technology for cement
        rscav(x,rgw,tp,sw)      Undiscovered resources available
        rsvav(x,rgw,tp,sw)      Proven reserves available
        rdflim(x,rgw,tp,sw)     Resource depletion limit
        prvlim(x,rgw,tp,sw)     Production-reserve limit
        carlev(rg,tp,sw)        Energy-related carbon emissions level - billion tons C
        ccemdf(rg,tp,sw)        Cement-related carbon emissions - billion tons C
        ch4lev(rg,tp,sw)        Energy-related methane emissions level - billion tons CH4
        othemitdef(ghg,rg,tp,sw) Non-energy Non-CO2 emissions definition
        ceqdef(rg,tp,sw)        Total emissions in billion tons C-eq
        wcardf(tp,sw)           Definition of world energy-related co2 emissions
        wch4df(tp,sw)           Definition of world energy-related ch4 emissions
        totemit(ghg,ctp,sw)     Total world emissions for each greenhouse gas
        clevbd(rg,tp,sw)        Upper bound on annual carbon emissions
        oggbd(ogg,rg,tp,sw)     Upper bound on annual ogg emissions
        ceqbd(rg,tp,sw)         Upper bound on total emissions in C-eq (using GWPs)
        itnbc(rg,sw)            Intertemporal budget constraint on net banked credits
        cbcdf(tp,rg,sw)         Definition of cumulative banked credits
*        affemitdef(tp,sw)
*        fracsum(sw)           sum of fractional aforestation programs
;

nweldf..    NWEL =E=  100 * 1000 * (

*        Utility for end-use model or consumption for aggregate model
           sum((pp,sw), prob(sw)* sum(rg, nwt(rg) * sum(swm(sw,sow,pp),udf(pp,rg)*
          ( log(U(rg,pp,sow))$tpclay(rg) + log(C(rg,pp,sow))$(not tpclay(rg))  ) )))

*        Damages for cost-benefit analysis
         + sum((ecp,sw), prob(sw)* sum(rg, nwt(rg) * sum(swm(sw,sow,ecp),udf(ecp,rg)*
            log(ELF(rg,ecp,sow))    ))) )

*        Penalty for RF (in special cases)
               - mu * RFMAX

*        Penalty for Temp
               - tau * TPMAX

*        Penalty for Conc.
               - zeta * CONMAX;

newdemand(tpclay(rg),st(pp,sw))..
        UN(rg,pp,sw) =e= uref(pp,rg) *
                (   alpha(pp,rg)  * (CN(rg,pp,sw)/cref(pp,rg))**rhot +
                 (1-alpha(pp,rg)) * (TN(rg,pp,sw)/tref(pp,rg))**rhot )**(1/rhot);

newprod(rg,st(pp,sw))..                YN(rg,pp,sw) =e= (
                                                         yref(pp,rg) *
        ((1-theta(pp,rg)) * ((KN(rg,pp,sw)/kref(pp,rg))**kpvs(rg) *
                (ln(pp,rg)/L(pp,rg))**(1-kpvs(rg)))**rho(pp,rg) +
            theta(pp,rg)  * ((EN(rg,pp,sw)/eref(pp,rg))**elvs(pp,rg) *
                (NN(rg,pp,sw)/nref(pp,rg))**(1-elvs(pp,rg)))**rho(pp,rg))**(1/rho(pp,rg))
                                                         )$%gross% +
                                                         (
                                                         ynref(pp,rg) *
        ((1-theta(pp,rg)) * ((KN(rg,pp,sw)/knref(pp,rg))**kpvs(rg) *
                             (ln(pp,rg)/ln(pp,rg))**(1-kpvs(rg)))**rho(pp,rg) +
            theta(pp,rg)  * ((EN(rg,pp,sw)/enref(pp,rg))**elvs(pp,rg) *
                             (NN(rg,pp,sw)/nnref(pp,rg))**(1-elvs(pp,rg)))**rho(pp,rg))**(1/rho(pp,rg))
                                                         )$%newvint%;
totalutl(tpclay(rg),st(pp(tp),sw))..
        UN(rg,tp,sw) =e= U(rg,tp,sw) - U(rg,tp-1,sw-sp(sw,tp))*speedt;

newcon(tpclay(rg),st(pp(tp),sw))..
        CN(rg,tp,sw) =e= C(rg,tp,sw) - C(rg,tp-1,sw-sp(sw,tp))*speedt;

newtrn(tpclay(rg),st(pp(tp),sw))..
        TN(rg,tp,sw) =e= T(rg,tp,sw) - T(rg,tp-1,sw-sp(sw,tp))*speedt;

newelec(rg,st(pp(tp),sw))..       EN(rg,tp,sw) =e= E(rg,tp,sw) - E(rg,tp-1,sw-sp(sw,pp))*speed(tp-1,rg);

newnon(rg,st(pp(tp),sw))..        NN(rg,tp,sw) =e= N(rg,tp,sw) - N(rg,tp-1,sw-sp(sw,pp))*speed(tp-1,rg);

* newcap(rg,st(pp,sw))..        KN(rg,pp,sw) =e= 5 * I(rg,pp,sw) + 5*speed(tp,rg)*I(rg,tp,sw-sp(sw,pp));
* New capital at start of period tp+1 should be a function of investment in period tp
newcap(rg,st(pp(tp),sw))..        KN(rg,tp,sw) =e= (((1 + grow(tp-1,rg))**10 - speed(tp-1,rg)) / (grow(tp-1,rg) + depr(rg)))
                                                    * I(rg,tp-1,sw-sp(sw,tp));

totalcap(rg,st(pp(tp),sw))..      KN(rg,tp,sw)  =e= K(rg,tp,sw)  - K(rg,tp-1,sw-sp(sw,tp))*speed(tp-1,rg);

totalprod(rg,st(pp(tp),sw))..     YN(rg,tp,sw) =e= Y(rg,tp,sw) -  Y(rg,tp-1,sw-sp(sw,tp))*speed(tp-1,rg);

termk(rg,tlast,sw)..            kterm(rg,sw) =e= K(rg,tlast,sw)*speed(tlast,rg)
*                                                + 5*speed(tlast,rg)*I(rg,tlast,sw);
                                  + (((1 + grow(tlast,rg))**10 - speed(tlast,rg)) / (grow(tlast,rg) + depr(rg)))
                                               * I(rg,tlast,sw);

*       Exogenous transportation demand

suptrans(trans(rg),pp,sw)..     sum(pvt, PT(pvt,pp,rg,sw)) =g= T(rg,pp,sw);

supliq(rg,pp,sw)$(not pfix(pp,rg))..      sum(lf, PR(lf,pp,rg,sw)) =g= PN("lqdu",pp,rg,sw)
                        + sum(v, htrt_v("oil-f",v,pp,rg) * PEV("oil-f",v,pp,rg,sw))
                        + sum(pvt, epvkt("lqtr",pvt,pp,rg) * PT(pvt,pp,rg,sw));

supvint(etv,rg,st(pp,sw))$(not pfix(pp,rg))..
                sum(vtp(v,pp), PEV(etv,v,pp,rg,sw)) =e= PE(etv,pp,rg,sw);

supelec(rg,st(pp,sw))..                sum(et, PE(et,pp,rg,sw)) =e=
         ( eveff(pp,rg) * PN("evhcl",pp,rg,sw) )$ntrns(rg)
        + (1/3.6) * sum(pvt, epvkt("elec",pvt,pp,rg) * PT(pvt,pp,rg,sw))$trans(rg)
        + E(rg,pp,sw);

supnon(rg,st(pp,sw))..  sum(nt, PN(nt,pp,rg,sw)) =e= N(rg,pp,sw);

concoal(rg,st(pp,sw))$(not pfix(pp,rg))..
        TPE("coal",pp,rg,sw)  =e=
*    Consumption of coal in electric sector
          sum(v, sum(xlc, htrt_v(xlc,v,pp,rg) * PEV(xlc,v,pp,rg,sw)))
*    Consumption of coal for liquid fuels
        + (1+syntpe) * PR("synf",pp,rg,sw)
*    Consumption of coal for direct use
        + PN("cldu",pp,rg,sw);

conoil(rg,st(pp,sw))$(not pfix(pp,rg))..
        TPE("oil",pp,rg,sw) =e= PR("petro",pp,rg,sw);

congas(rg,st(pp,sw))$(not pfix(pp,rg))..
        TPE("gas",pp,rg,sw)  =e=
*    Consumption of gas in electric sector
          sum(v, sum(xlg, htrt_v(xlg,v,pp,rg) * PEV(xlg,v,pp,rg,sw)))
*    Consumption of gas in passenger vehicles
        + sum(pvt, epvkt("cng",pvt,pp,rg) * PT(pvt,pp,rg,sw))$trans(rg)
*    Consumption of gas for direct use
        + PN("gsdu",pp,rg,sw);

connuc(rg,st(pp,sw))$(not pfix(pp,rg))..
         TPE("ur",pp,rg,sw) =e=  10 * PE("nuc-1",pp,rg,sw);

supx(x,rg,st(pp,sw))$xmap(x,rg)..
        EX(x,pp,rg,sw) - sum(xtrd(x,trd), NTX(trd,pp,rg,sw)) =g= TPE(x,pp,rg,sw);

wsupx(x,st(pp,sw))$xmap(x,"world")..
        EX(x,pp,"world",sw) =g= sum(rg, TPE(x,pp,rg,sw));

cxdef(xmap(x,rgw),st(tp,sw))$ursc0(x,"2000",rgw)..
         CX(x,tp,rgw,sw)  =e= 1 - (URSC(x,tp,rgw,sw)+PRSV(x,tp,rgw,sw))/(ursc0(x,"2000",rgw) + prsv0(x,"2000",rgw));

xcstdef(xmap(x,rgw),st(pp,sw))..
         XCST(x,pp,rgw,sw) =e= xcst0(x,rgw) + (xcst_u(x)-xcst0(x,rgw))*CX(x,pp,rgw,sw);

cadef(xmap(x,rgw),st(tp,sw))$ursc0(x,"2000",rgw)..
         CA(x,tp,rgw,sw)  =e= 1 - URSC(x,tp,rgw,sw) /ursc0(x,"2000",rgw);

acstdef(xmap(x,rgw),st(pp,sw))..
         ACST(x,pp,rgw,sw) =e= acst0(x,rgw) + (acst_u(x)-acst0(x,rgw))*CA(x,pp,rgw,sw);

cgdef(st(tp+1,sw))..
         CG(tp+1,sw) =e= CG(tp,sw-sp(sw,tp+1)) + 5 * sum(rg, TPE("coal",tp,rg,sw-sp(sw,tp+1)) +
                                                             TPE("coal",tp+1,rg,sw) );

scdef(rg,st(tp+1,sw))..
         SC(tp+1,rg,sw) =e= SC(tp,rg,sw-sp(sw,tp+1)) +
                  5 * (sum((v,ccs), (cprt(ccs,v,rg)/(1-cprt(ccs,v,rg))) * cece_v(ccs,v,tp+1,rg) * PEV(ccs,v,tp+1,rg,sw) +
                                    (cprt(ccs,v,rg)/(1-cprt(ccs,v,rg))) * cece_v(ccs,v,tp,rg) * PEV(ccs,v,tp,rg,sw-sp(sw,tp+1))) +
                       cec_cem * 0.9 * (CEMCCS(tp+1,rg,sw) + CEMCCS(tp,rg,sw-sp(sw,tp+1))));
*----------------------------------------------------------------------------------------------------------------------------
*       Natural gas is limited to supplying 50% of the electric energy market.

egfrac(rg,st(pp,sw))$(not pfix(pp,rg))..   sum(xlg(et), PE(et,pp,rg,sw)) =l= 0.5 * E(rg,pp,sw);

*       Coal is limited to supplying 60% of the electric energy market in US (other fraction elsewhere).

ecfrac(rg,st(pp,sw))$(not pfix(pp,rg))..   sum(xlc(et), PE(et,pp,rg,sw)) =l= coallim(rg) * E(rg,pp,sw);

*   CCS is limited to supplying 100% of the electric energy market.

ccsfrac(rg,st(pp,sw))$(not pfix(pp,rg))..   sum(ccs(et), PE(et,pp,rg,sw)) =l=  1 * E(rg,pp,sw);

*       Wind and solar-lc (and existing intermittents) are collectively limited to 20% of the electric market:

intfrac(rg,st(pp,sw))$(not pfix(pp,rg))..         sum(irnw(et), PE(et,pp,rg,sw)) =l= 0.2 * E(rg,pp,sw);

*       Constrained (i.e. non-backstop) renewables are limited to 30% of electric market:

constrnw(rg,st(pp,sw))$(not pfix(pp,rg))..        sum(crnw(et), PE(et,pp,rg,sw)) =l= 0.3 * E(rg,pp,sw);

*   Total non-hydro non-waste renewable production has an exogenous limit in the US:

rnwusa(st(pp,sw))..     sum(usa(rg), sum(rnw(et), PE(et,pp,rg,sw))) =l= rnwusacap(pp);

*   Include RPS in some scenarios
*-----------------------------------------------------------------------------------------------------------------------
$if %rps%==yes $include "RPS.tab"

*   Nuclear power in USA in limited pessimistic case does not exceed existing production level:

nucusa(st(pp,sw))..     sum(usa(rg), sum(nuc, PE(nuc,pp,rg,sw))) =l= nucusacap(pp);

*  ETAC limit on 2010-20 nuclear builds in US

nucusa20(sw)..     sum(usa(rg), PEV("nuc-1","2020","2020",rg,sw)) =l= 0.189;

*       Natural gas is limited to 50% of the non-electric energy market.

gfrac(rg,st(pp,sw))$(not pfix(pp,rg))..    PN("gsdu",pp,rg,sw) =l= .5 * N(rg,pp,sw);

*       Biofuels and synfuels are limited to 50% of the liquids refining market.

* lffrac(nlf,rg,st(pp,sw))..      PR(nlf,pp,rg,sw) =l= 0.5 * sum(llf, PR(llf,pp,rg,sw));


costtrn(trans(rg),st(tp,sw))..          TC(rg,tp,sw) =e= 0.001 * (

*       Non-energy cost of transportation services

          sum(pvt, vcst(pvt,tp,rg) * PT(pvt,tp,rg,sw))

*       Wholesale energy costs for oil, gas, and electricity are captured upstream
*       lcst() and tfcst() are only non-zero for fuels without markets elsewhere
*       i.e. synfuels, biofuels, and non-elec backstop

        + sum(tf, tfcst(tf,tp,rg)*sum(pvt, epvkt(tf,pvt,tp,rg)*PT(pvt,tp,rg,sw)))

*       For end-use transportation service, must include retail margins and taxes

        + sum(tf, retail(tf,tp,rg)*sum(pvt, epvkt(tf,pvt,tp,rg)*PT(pvt,tp,rg,sw))) );

*   The final equation in the ETA submodel defines EC, the aggregate costs
*   incurred within the energy sector.  This is one of the competing claims on
*   aggregate output.

costnrg(rg,st(tp,sw))$pp(tp)..    EC(rg,tp,sw) =g= 0.001 * (

*   Costs for electric generation:

        sum((etv, vtp(v,tp)), PEV(etv,v,tp,rg,sw) * ecst_v(etv,v,tp,rg))

      + sum(et$(not etv(et)), PE(et,tp,rg,sw)*ecst(et,tp,rg))

*    Costs for oil and gas production (avoid double counting in one-world model)

      + sum(x,  EX(x,tp,rg,sw) * XCST(x,tp,rg,sw))$(not sameas(rg,"world"))

*    Costs for oil and gas exploration

      + sum(x,  RA(x,tp,rg,sw) * ACST(x,tp,rg,sw))$(not sameas(rg,"world"))

*    Production and exploration costs are shared among regions according to consumption for world supply model:

      + sum(x$xmap(x,"world"), (EX(x,tp,"world",sw) * XCST(x,tp,"world",sw)
                              + RA(x,tp,"world",sw) * ACST(x,tp,"world",sw)) * xshr(x,tp,rg,sw))

*    Costs for petroluem refining (incurred by consuming region)

      + PR("petro",tp,rg,sw) * rcst(tp,rg)

*    Costs for production plus refining of SYNF and BFUEL:

      + sum(nlf(lf), PR(lf,tp,rg,sw)*lcst(lf,tp,rg))

*    Costs for other non-electric technologies

      + sum(nt, PN(nt,tp,rg,sw)*ncst(nt,tp,rg))

*   Allowance for oil-gas (intraregion pipeline transport) price differential

                    + ogpd(rg)*TPE("gas",tp,rg,sw)

*   Taxes on electricity, nonelectric energy, carbon, and methane.

                    + xntax(rg)*N(rg,tp,sw)

                    + xetax(rg)*E(rg,tp,sw)

                    + ghgtax("co2",tp,rg,sw)*CLEV(tp,rg,sw)

                    + ghgtax("ch4",tp,rg,sw)*MLEV(tp,rg,sw)

                    + rsktax("n",tp,rg,sw)*sqr(sum(nuc, PE(nuc,tp,rg,sw)))

                    + rsktax("c",tp,rg,sw)*sqr(sum(ccs, PE(ccs,tp,rg,sw)))

                    + rsktax("p",tp,rg,sw)*sqr(PE("nuc-adv",tp,rg,sw))

*   Credit for lump-sum rebate of tax revenues.

                    - taxrev(tp,rg,sw)

*   Inecurity price premium on fuel exports

                    + sum(trd, ipprem(trd,tp,rg)*EXPRT(trd,tp,rg,sw))

*   Credit for insecurity premium revenue recycling

                    - iprev(tp,rg,sw)


*   Transportation costs for interregional trade.

                    + sum(trd, cstexp(trd)*EXPRT(trd,tp,rg,sw)) );

* Cost of non-energy related abatement:

costabt(rg,st(pp,sw))..    AC(rg,pp,sw) =g= 0.001 * (
*   Cost of providing abatement of non-energy related OGG.
                    + sum ((ghg,abx), 0.001 * abcst(abx)*ABATE(ghg,abx,pp,rg,sw))
*   Cost of using capture technology for cement production
                    + CEMCCS(pp,rg,sw) * 0.9 * cec_cem * cemccs_cost(pp,rg)
*   Cost of abating bunker fuel emissions:
                    + sh(pp,rg) * OTHCABT(pp,sw) * bunker_abcst(pp)
*   Regional cost of global afforestation - based on share parameter
*                    + sh(pp,rg) *  acf * SUM((sy,ps), price(sy,pp,ps)*cab(sy,pp,ps)*AFF(sy,ps,sw))
);

* Payments of non-energy related emissions taxes (or subsidy for afforestation)
* For feasibility, payments are only made on emissions that can be abated

nontaxdef(rg,st(pp,sw))..

         NONTAX(rg,pp,sw) =e= 0.001 * (
                 sum(ogg(ghg), ghgtax(ghg,pp,rg,sw) * (OTHEMIT(ghg,pp,rg,sw) - sum(ogs, oggmin(ghg,ogs,pp,rg))))
                 + ghgtax("co2",pp,rg,sw) * CCEM(pp,rg,sw)
*                 + ghgtax("co2",pp,rg,sw) * sh(pp,rg) * AFFEMIT(pp,sw)
);

*   Competing claims on resources include consumption, investment, energy costs,
*   transportation costs, non-CO2 abatement costs, market damages and net exports of the composite numeraire good.

cc(rg,st(pp,sw))..
        Y(rg,pp,sw) =e= C(rg,pp,sw) + I(rg,pp,sw) + EC(rg,pp,sw) + TC(rg,pp,sw) + AC(rg,pp,sw) + NONTAX(rg,pp,sw) + MD(rg,pp,sw) + NTX("nmr",pp,rg,sw);
*
*   Net exports of interregionally tradeable goods must be balanced with net
*   imports. (Not true in base year for oil and gas)

trdbal(trd,st(pp,sw))..  sum(rg, NTX(trd, pp, rg,sw)) =e= 0;

*   In accounting for total energy costs, interregional transport costs are
*   proportional to net exports. It will be optimal for EXPRT variables to be
*   zero when NTX < 0.

xptdef(trd,rg,st(tp,sw))..   NTX(trd,tp,rg,sw) =l= EXPRT(trd,tp,rg,sw);

iptdef(trd,rg,st(tp,sw))..   IMPRT(trd,tp,rg,sw) =e= EXPRT(trd,tp,rg,sw) - NTX(trd,tp,rg,sw);

*   Pipeline exports must be received by connected countries
pplxpt(rg,st(pp,sw)).. EXPRT("ngpl",pp,rg,sw) =l= sum(ppl(pp,rrg,rg), IMPRT("ngpl",pp,rrg,sw));

*   Pipeline imports must be sent from connected countries
pplipt(rg,st(pp,sw)).. IMPRT("ngpl",pp,rg,sw) =l= sum(ppl(pp,rg,rrg), EXPRT("ngpl",pp,rrg,sw));
*-------------------------------------------------------------------------------------------------------------------
*   New vintage electric technologies are subject to decline limits.

decev(rg,dle(et),vnw(v),st(tp+1,sw))$(vtp(v,tp) and not pfix(tp+1,rg))..

        PEV(et,v,tp+1,rg,sw) =g= PEV(et,v,tp,rg,sw-sp(sw,tp+1))*(decf(tp,rg)**10);

*   Existing post-90 vintages of fossil generation must decline.

fosdec(etf,v,rg,st(pp(tp),sw))$(dlv(v) and vex(v) and vtp(v,tp-1) and not pfix(tp,rg))..

        PEV(etf,v,tp,rg,sw) =l= PEV(etf,v,tp-1,rg,sw-sp(sw,tp))*(decf(tp-1,rg)**10);

*   Nuclear power supply must be monotonically increasing until retirement (both gen3 and gen4, separately)

nucdec(nuc,dlv(v),rg,st(tp+1,sw))$(vtp(v,tp) and not pfix(tp+1,rg))..
        PEV(nuc,v,tp+1,rg,sw) =g= PEV(nuc,v,tp,rg,sw-sp(sw,tp+1))$vnuc(v,tp+1);

dece(rg,dle,st(tp+1,sw))$(not etv(dle) and not pfix(tp+1,rg))..
        PE(dle,tp+1,rg,sw) =g= PE(dle,tp,rg,sw-sp(sw,tp+1))*decf(tp,rg)**10;

decn(rg,dln,st(pp(tp),sw))$(not pfix(tp,rg))..
        PN(dln,tp,rg,sw) =g= PN(dln,tp-1,rg,sw-sp(sw,tp))*decf(tp-1,rg)**10;

dect(trans(rg),pvt,st(tp+1,sw))$(not pfix(tp+1,rg))..
        PT(pvt,tp+1,rg,sw) =g= PT(pvt,tp,rg,sw-sp(sw,tp+1)) * decf(tp,rg)**10;

decr(rg,nlf(lf),st(tp+1,sw))$(not pfix(tp+1,rg))..
        PR(lf,tp+1,rg,sw) =g= PR(lf,tp,rg,sw-sp(sw,tp+1)) * decf(tp,rg)**10;

decx(rg,x,st(tp+1,sw))$xmap(x,rg)..
        EX(x,tp+1,rg,sw) =g= EX(x,tp,rg,sw-sp(sw,tp+1)) * 0.9**10;
*--------------------------------------------------------------------------------------------------------------
*   There are limits on the expansion rates of specific technologies.

*   Vintages cannot expand in subsequent years
expev(rg,et,xlv(v),st(tp+1,sw))$(vtp(v,tp) and not pfix(tp+1,rg))..
        PEV(et,v,tp+1,rg,sw) =l= PEV(et,v,tp,rg,sw-sp(sw,tp+1));

*   New vintage Nuclear expansion in Annex B (refers to next generation)

expnucanb(nuc,anb(rg),st(tp+1,sw))$(not pfix(tp+1,rg))..
        sum(vtp(vnw,tp+1), PEV(nuc,vnw,tp+1,rg,sw)) =l= 3 * sum(vtp(vnw,tp), PEV(nuc,vnw,tp,rg,sw-sp(sw,tp+1)))
        + 0.08 * E(rg,tp+1,sw);

*  Expansion of nuclear in non-Annex B applies to all capacity

expnucnanb(nuc,nanb(rg),st(tp+1,sw))$(not pfix(tp+1,rg))..
        PE(nuc,tp+1,rg,sw) =l= expf(rg)**10 * PE(nuc,tp,rg,sw-sp(sw,tp+1)) + nshf(rg) * E(rg,tp+1,sw);

*   Expansion of coal- and gas-fired electricity (applies to all capacity, not just existing):

expc(rg,st(tp+1,sw))$(not pfix(tp+1,rg))..
        sum(xlc, PE(xlc,tp+1,rg,sw)) - nshf(rg) * E(rg,tp+1,sw)
                =l=   expf(rg)**10 * sum(xlc, PE(xlc,tp,rg,sw-sp(sw,tp+1)));

expg(rg,st(tp+1,sw))$(not pfix(tp+1,rg))..
        sum(xlg, PE(xlg,tp+1,rg,sw)) - nshf(rg) * E(rg,tp+1,sw)
                =l=   expf(rg)**10 * sum(xlg, PE(xlg,tp,rg,sw-sp(sw,tp+1)));

* Expansion of renewables (applies to each technology individually, including existing capacity)

exprnw(rnw,rg,st(tp+1,sw))$(not pfix(tp+1,rg))..
        PE(rnw,tp+1,rg,sw) =l= expf(rg)**10 * PE(rnw,tp,rg,sw-sp(sw,tp+1))
        + nshf(rg) * E(rg,tp+1,sw);

* Expansion of CCS applies to sum of coal and gas, follows faster nuclear rate

expccs(rg,st(tp+1,sw))$(not pfix(tp+1,rg))..
        sum(ccs(et),PE(et,tp+1,rg,sw)) =l= 3 * sum(ccs(et),PE(et,tp,rg,sw-sp(sw,tp+1)))
                + 0.08 * E(rg,tp+1,sw);
*------------------------------------------------------------------------------------------------------------------------
*       Expansion of non-electric energy:

expn(rg,xln,st(tp+1,sw))$(not pfix(tp+1,rg))..
        PN(xln,tp+1,rg,sw) =l= nxpf(rg)**10 * PN(xln,tp,rg,sw-sp(sw,tp+1))
        + nshf(rg) * N(rg,tp+1,sw) ;


evfrac(ntrns(rg),st(pp,sw))$(not pfix(pp,rg))..
        PN("evhcl",pp,rg,sw) =l= evlimit(pp) * N(rg,pp,sw);

expt(trans(rg),xlt(pvt),st(tp+1,sw))$(not pfix(tp+1,rg))..
        PT(pvt,tp+1,rg,sw) =l= expf(rg)**10 * PT(pvt,tp,rg,sw-sp(sw,tp+1))
        + nshf(rg) * T(rg,tp+1,sw);

expr(rg,nlf(lf),st(tp+1,sw))$(not pfix(tp+1,rg))..
        PR(lf,tp+1,rg,sw) =l= expf(rg)**10 * PR(lf,tp,rg,sw-sp(sw,tp+1))
        + nshf(rg) * sum(llf, PR(llf,tp+1,rg,sw));

expcm(rg,st(tp+1,sw))$(not pfix(tp+1,rg))..
        CEMCCS(tp+1,rg,sw) =l= expf(rg)**10 * CEMCCS(tp,rg,sw) + nshf(rg) * cement(tp+1,rg);

* Only apply LNG expansion constraints when pipeline option is turned on
$if not set ngpl $goto resources

equations
         explngx(rg,tp,sw)       Expansion rate of LNG exports
         explngm(rg,tp,sw)       Expansion rate of LNG exports;

explngx(rg,st(tp+1,sw))..
        EXPRT("lng",tp+1,rg,sw) =l= lngx(rg)**10 * EXPRT("lng",tp,rg,sw-sp(sw,tp+1))
        + lngf(rg) * TPE("gas",tp+1,rg,sw);

explngm(rg,st(tp+1,sw))..
        IMPRT("lng",tp+1,rg,sw) =l= lngx(rg)**10 * IMPRT("lng",tp,rg,sw-sp(sw,tp+1))
        + lngf(rg) * TPE("gas",tp+1,rg,sw);

$label resources

*   Changes in remaining undiscovered exhaustible resources are determined by a
*   distributed lag function of reserve additions.

rscav(xmap(x,rgw),st(tp+1,sw))..
        URSC(x,tp+1,rgw,sw) + .5*10 * RA(x,tp+1,rgw,sw) =e=
                URSC(x,tp,rgw,sw-sp(sw,tp+1)) -.5*10 * RA(x,tp,rgw,sw-sp(sw,tp+1));

*   Changes in remaining proven reserves are determined by a distributed lag
*   function of reserve additions less production of exhaustible resources.

rsvav(xmap(x,rgw),st(tp+1,sw))..
        PRSV(x,tp+1,rgw,sw) - .5*10 * (RA(x,tp+1,rgw,sw) - EX(x,tp+1,rgw,sw)) =e=
        PRSV(x,tp,rgw,sw-sp(sw,tp+1))
                + .5*10 * (RA(x,tp,rgw,sw-sp(sw,tp+1)) - EX(x,tp,rgw,sw-sp(sw,tp+1)));

*   Reserve additions cannot exceed a fixed fraction of remaining resources.

rdflim(xmap(x,rgw),st(pp,sw))..
        RA(x,pp,rgw,sw)    =l= rdf(x,rgw) * URSC(x,pp,rgw,sw);

*   Production cannot exceed a fixed fraction of proven reserves.

prvlim(xmap(x,rgw),st(pp,sw))..
         EX(x,pp,rgw,sw) =l= prv(x,rgw) * PRSV(x,pp,rgw,sw);

$ontext

Carbon emissions are determined by the carbon coefficients and the production
levels for each electric and nonelectric technology. Emissions may
be reduced by high-cost carbon limit relaxation activities.

$offtext

carlev(rg,st(pp,sw))..
        CLEV(pp,rg,sw) =e=
              sum((etv,vtp(v,pp)), cece_v(etv,v,pp,rg) * PEV(etv,v,pp,rg,sw))
            + sum(lf, cecl(lf,pp) * PR(lf,pp,rg,sw))
            + cec("coal") * PN("cldu",pp,rg,sw)
            + cec("gas")  * sum(pvt, epvkt("cng",pvt,pp,rg) * PT(pvt,pp,rg,sw))$trans(rg)
            + cec("gas")  * PN("gsdu",pp,rg,sw)

*       Adjustment for base year inconsistency of IEA energy data and ORNL emissions data:

            + carbdif(rg);

ccemdf(rg,st(pp,sw))..
         CCEM(pp,rg,sw) =e= cec_cem * (cement(pp,rg) - 0.9 * CEMCCS(pp,rg,sw));

ch4lev(rg,st(pp,sw))..
         MLEV(pp,rg,sw) =e= 0.001 * sum(fos, ch4rate(fos,pp,rg)*TPE(fos,pp,rg,sw));

* Non-energy-related emissions of other GHGs (convert from carbon equivalent)

othemitdef(ogg(ghg),rg,tp,sw)..
        OTHEMIT(ghg,tp,rg,sw) =e= (sum(ogs, oggbline(ghg,ogs,tp,rg)) -
                                   sum(abx, ABATE(ghg,abx,tp,rg,sw)))/(1000 *gwp(ghg));

* Total emissions by region expressed in C-eq terms

ceqdef(rg,tp,sw)..
        CEQ(tp,rg,sw) =e= CLEV(tp,rg,sw) + CCEM(tp,rg,sw) + MLEV(tp,rg,sw)*gwp("ch4") +
                          sum(ogg, OTHEMIT(ogg,tp,rg,sw)*gwp(ogg));

* Global total energy-related carbon emissions:

wcardf(pp(tp),sw)$st(pp,sw)..   EM("co2",tp,sw) =e= sum(rg, CLEV(tp,rg,sw))
*                               includes international bunkers (less abatement)
                                + bunker(tp) - OTHCABT(tp,sw)
*                               and adjustment to match global ORNL total
                                + gcarbdif;

* Global total energy-related methane emissions

wch4df(pp(tp),sw)$st(pp,sw)..   EM("ch4",tp,sw) =e= sum(rg, MLEV(tp,rg,sw));

$ontext
* Global total net land emissions

landemit(st(tp,sw))..          NETLAND(tp,sw) =e= luc(tp) + CO2FERT(tp);

* Global total terrestrial sink from CO2 fertilization term (from MAGICC)

terrsink(st(tp,sw))..           CO2FERT(tp) =e=
$offtext
* Global total emissions by GHG

totemit(ghg,tp,sw)..
      TOTEM(ghg,tp,sw)  =e=     EM(ghg,tp,sw)
                                + sum(rg, CCEM(tp,rg,sw))$sameas(ghg,"co2")
                                + NETLAND(tp,sw)$sameas(ghg,"co2")
*                               + AFFEMIT(tp,sw)$sameas(ghg,"co2")
                                + sum(rg, OTHEMIT(ghg,tp,rg,sw));

* Exogenous constraint on carbon emissions:

clevbd(rg,st(pp,sw))$(not pfix(pp,rg))..
         CLEV(pp,rg,sw) + CCEM(pp,rg,sw) =L=  carlim(pp,rg) - NTX("crt-co2",pp,rg,sw) - NBC(pp,rg,sw);

* Exogenous constraint on other GHGs:

oggbd(ogg,rg,st(pp,sw))$(not pfix(pp,rg))..
        ghgtol(ogg) * (OTHEMIT(ogg,pp,rg,sw) + MLEV(pp,rg,sw)$sameas(ogg,"ch4")) =L=
        ghgtol(ogg) * (ogglim(ogg,pp,rg) - sum(xghg(ogg,trd),NTX(trd,pp,rg,sw)));

* In some cases, quantity constraints are enforced in C-eq terms:

ceqbd(rg,st(pp,sw))$(not pfix(pp,rg))..
         CEQ(pp,rg,sw) =L= ceqlim(pp,rg) - NTX("crt-ceq",pp,rg,sw) - NBC(pp,rg,sw);


* Intertemporal budget constraint on banked credits

itnbc(rg,sw).. sum(tnbc(tp,rg), NBC(tp,rg,sw)) =e= 0;

* Cumulative banked credits - must be positive if borrowing is not allowed:

cbcdf(tp,rg,sw).. CBC(tp,rg,sw) =e= sum(tnbc(ttp,rg)$(yr(ttp) le yr(tp)), NBC(ttp,rg,sw));

* * * * * * * * * * * * * Bounds and initialization * * * * * * * * * * * * *

*       Assign default initial values for variables which enter nonlinear
*       functions in the model (level values for linear variables are not
*       required):

K.L(rg,pp,sw)    =  kref(pp,rg);
Y.L(rg,pp,sw)    =  yref(pp,rg);
C.L(rg,tp,sw)    =  cref(tp,rg);
I.L(rg,pp,sw)    =  sum(tbase, i0(tbase,rg)*l(pp,rg));

CN.L(tpclay(rg),pp,sw) = C.L(rg,pp,sw) * ln(pp,rg)/l(pp,rg);
TN.L(tpclay(rg),pp,sw) = tref(pp,rg) * ln(pp,rg)/l(pp,rg);
U.L(tpclay(rg),tp,sw) = C.L(rg,tp,sw)/alpha(tp,rg);
UN.L(tpclay(rg),pp,sw) = U.L(rg,pp,sw)* ln(pp,rg)/l(pp,rg);

KN.L(rg,pp,sw)   = sum(tbase, k0(tbase,rg) * ln(pp,rg)/l(pp,rg));
YN.L(rg,pp,sw)   = sum(tbase, y0(tbase,rg) * ln(pp,rg)/l(pp,rg));
EN.L(rg,pp,sw)   = sum(tbase, e0(tbase,rg) * ln(pp,rg)/l(pp,rg));
NN.L(rg,pp,sw)   = sum(tbase, n0(tbase,rg) * ln(pp,rg)/l(pp,rg));
E.L(rg,pp,sw)   =  eref(pp,rg);

TPE.L("ur",pp,rg,sw) = 10 * pe0("nuc-1","2000",rg);
rsktax(rsk,pp,rg,sw) = (wtp0(rsk,rg)/2) * (1/ (refshr(rsk,rg)*E.L(rg,pp,sw)));
MD.L(rg,ctp,sw) = 0;
ELF.L(rg,ctp,sw) = 1;


NWEL.L = 0;
RFMAX.L = 0;
TPMAX.L = 0;
CLEV.L(tp,rg,sw) = 0;

* Unless ptefx option is selected, climate model is only evaluated
* over first century in initial iteration:

cpp(ctp)$(not ctfix(ctp)) = yes;
$if not set ptefx cpp(ctp)$(not tp(ctp)) = no;
* Define economic-climate projection periods
ecp(cpp) = yes; loop(tfix, ecp(ctp)$sameas(ctp,tfix) = no;);

*       Read the climate model and any additional equations
*       which may be part of the previous basis:
*-------------------------------------------------------------------------------------------------
*QUI SI INIZIA A FAR SCENARI 
$if set contrg          $set climate yes
$if set rfktrg          $set climate yes
$if set tmptrg          $set climate yes
$if set inctrg          $set climate yes
$if set tgctrg          $set climate yes
$if %mktdam%==yes       $set climate yes
$if %nmktdam%==yes      $set climate yes
$if set mu              $set climate yes
$if set tau             $set climate yes
$if set zeta            $set climate yes
$if set climate_fb      $set climate yes
$if set climate $include "climate.gms"

*       Read an advanced basis if one is specified.  A basis includes
*       the point values as well as the Negishi weights:

$if exist %basdir%%basis%.gdx $goto readbasis
$if exist %basdir%m7_p.gdx $goto defaultbasis

*       No basis is read -- need to have an initial value:

$if set climate TOTEM.L(ghg,ctp,sw) = 1;

$goto bounds

* These parameters are unloaded and need to be read from %run%_nwt.gdx:
* execute_unload 'm7_nwt',nwt, taxrev, rsktax, ghgsbd, nashsbd, iprev, kterm;

$label readbasis
execute_loadpoint '%basdir%%basis%.gdx';
execute_load      '%basdir%%basis%_nwt.gdx',nwt,taxrev,rsktax,ghgsbd,iprev,kterm;
$if set nash execute_load '%basdir%%basis%_nwt.gdx',nashsbd; 
* $if set nash parameter nashsbd_(tp,rg,sw); execute_load '%basdir%%basis%_nwt.gdx',nashsbd_=nashsbd; nashsbd("co2",tp,rg,sw) = nashsbd_(tp,rg,sw);
execute_load      '%basdir%%basis%.gdx',TOTEM.L;
$if set ptefx TOTEM.FX(ghg,ctp,sw)$(not tp(ctp)) = TOTEM.L(ghg,ctp,sw);

$goto bounds

$label defaultbasis
execute_loadpoint '%basdir%m7_p.gdx';
execute_load      '%basdir%m7_nwt.gdx',nwt,taxrev,rsktax,ghgsbd,nashsbd,iprev,kterm;

$label bounds

*       Fix post-terminal emissions:  feasibility sensitive to this choice!
*       The fixed pathway will be subsequently refined in the iterative loop:

* $if set climate TOTEM.FX(ghg,ctp,sw)$(not tp(ctp)) = TOTEM.L(ghg,ctp,sw);

*  This section of the code contains the base year benchmarks of
*  variables and other variable bounds.

E.FX(rg,tfix,sw) = e0(tfix,rg);
N.FX(rg,tfix,sw) = n0(tfix,rg);
Y.FX(rg,tfix,sw) = y0(tfix,rg);
K.FX(rg,tfix,sw) = k0(tfix,rg);
C.FX(rg,tfix,sw) = c0(tfix,rg);
I.FX(rg,tfix,sw) = i0(tfix,rg);

U.FX(rg,tfix,sw) = U.L(rg,tfix,sw);
T.FX(rg,tfix,sw) = tref(tfix,rg);

* The following bounds help to avoid nasty program calls.

parameter        lotol /0.05/;

K.LO(rg,pp,sw)   =  lotol*kref(pp,rg);

Y.LO(rg,pp,sw)   =  lotol*yref(pp,rg);

C.LO(rg,pp,sw)   =  max(0.1, lotol*C.L(rg,pp,sw));
I.LO(rg,pp,sw)   =  lotol*I.L(rg,pp,sw);
E.LO(rg,pp,sw)   =  lotol*eref(pp,rg);
N.LO(rg,pp,sw)   =  lotol*nref(pp,rg);

KN.LO(rg,pp,sw)  =  lotol*KN.L(rg,pp,sw);
YN.LO(rg,pp,sw)  =  lotol*YN.L(rg,pp,sw);
EN.LO(rg,pp,sw)  =  lotol*EN.L(rg,pp,sw);
NN.LO(rg,pp,sw)  =  lotol*NN.L(rg,pp,sw);

CN.LO(tpclay(rg),pp,sw) = lotol*CN.L(rg,pp,sw);
TN.LO(tpclay(rg),pp,sw) = lotol*TN.L(rg,pp,sw);
U.LO(tpclay(rg),pp,sw) = lotol*U.L(rg,pp,sw);
UN.LO(tpclay(rg),pp,sw) = lotol*UN.L(rg,pp,sw);

* Electric energy initialization and bounds

PE.FX(et,tfix,rg,sw) =  pe0(et,tfix,rg);
PEV.FX(etv,v,tfix,rg,sw) = pev0(etv,v,tfix,rg);

* ECAP applies only to new technologies (no deployment today)
PE.UP(etn,pp,rg,sw) = ecap(etn,pp,rg);
* Specific limit for wind
PE.UP("wind",pp,rg,sw) = sum(xrg(rg,sub), wind_s(sub,pp));

* Hydroelectric is fixed to WEO levels through 2030

PE.FX("hydro",pp,rg,sw)  =  weo("hydro",pp,rg);
PE.FX("hydro",pp,rg,sw)$(yr(pp) > 2030) = weo("hydro","2030",rg);


* Other renewables remain fixed at projected levels through 2010

PE.FX("oth-r",pp,rg,sw) = pe0("oth-r",pp,rg);
PE.FX("oth-r",pp,rg,sw)$(yr(pp) > 2010) = pe0("oth-r","2010",rg);


* Total nuclear is fixed to WEO levels till 2010 for OECD, till 2030 outside

PE.FX("nuc-1",pp,oecd,sw)$(yr(pp) <= 2010)  = weo("nuc-1",pp,oecd);
PE.FX("nuc-1",pp,noecd,sw)$(yr(pp) <= 2030) = weo("nuc-1",pp,noecd);

* For pessimistic nuclear scenarios, fix new vintage nuc-1 to zero:

$if %nuc%==pess PEV.UP("nuc-1",v,pp,oecd,sw)$(vyr(v) > 2010) = 0;
$if %nuc%==pess PEV.UP("nuc-1",v,pp,noecd,sw)$(vyr(v) > 2030) = 0;

*  Fix retired nuclear vintages to zero:

PEV.FX(nuc,v,tp,rg,sw)$(not vnuc(v,tp)) = 0;

* Oil generation observed levels for 2000, 2010 are treated as
* upper bounds for OECD, fixed points for Non-OECD:

PE.UP("oil-f",pp,oecd,sw)$(yr(pp) <=2010)  = pe0("oil-f",pp,oecd);
PE.FX("oil-f",pp,noecd,sw)$(yr(pp) <=2010) = pe0("oil-f",pp,noecd);

* New vintages of oil are not allowed after 2010:

PEV.UP("oil-f",v,pp,rg,sw)$(vyr(v) > 2010) = 0;

* Assume all oil generation is phased out immediately after 2010 in OECD:
* Non-OECD oil generation declines as a dle technology

PE.UP("oil-f",pp,oecd,sw)$(yr(pp) > 2010) = 0;

* 1990- vintage capacity must be retired according to a fixed schedule
* For some technologies / regions, early retirement is allowed

loop(tbase,
* No early retirement of nuclear
         PEV.FX("nuc-1","1990",pp,rg,sw) = pev0("nuc-1","1990",tbase,rg) * ertr90("nuc-1",pp,rg);
*----------------------------------------------------------
*QUI 
* No early retirement of oil (goes to zero after 2010 in OECD)
         PEV.FX("oil-f","1990",pp,rg,sw) = pev0("oil-f","1990",pp,rg) * ertr90("oil-f",pp,rg);
* No early retirement of gas and coal outside OECD
         PEV.FX("gas-f","1990",pp,noecd,sw) = pev0("gas-f","1990",tbase,noecd) * ertr90("gas-f",pp,noecd);
         PEV.FX("coal-f","1990",pp,noecd,sw) = pev0("coal-f","1990",tbase,noecd) * ertr90("coal-f",pp,noecd);
* Early retirement of gas and coal allowed inside OECD
         PEV.UP("gas-f","1990",pp,oecd,sw) = pev0("gas-f","1990",tbase,oecd) * ertr90("gas-f",pp,oecd);
         PEV.UP("coal-f","1990",pp,oecd,sw) = pev0("coal-f","1990",tbase,oecd) * ertr90("coal-f",pp,oecd);
* No retirement of wind (ertr90 = 1)
         PEV.FX("wind","1990",pp,rg,sw) = pev0("wind","1990",tbase,rg) * ertr90("wind",pp,rg);
);



*  Fix inactive vintages to zero:

PEV.FX(et,v,tp,rg,sw)$(not (vtp(v,tp) and etv(et))) = 0;


* Non-electric energy initialization and bounds

PN.FX(nt,tfix,rg,sw) = pn0(nt,tfix,rg);
PR.FX(lf,tfix,rg,sw) = pr0(lf,tfix,rg);

PN.UP(nt,pp,rg,sw) = ncap(nt,pp,rg);
PR.UP(lf,pp,rg,sw) = rcap(lf,pp,rg);

TPE.FX(fos,tfix,rg,sw) = foscon0(fos,tfix,rg);
TPE.FX("ur",tfix,rg,sw) = 10 * pe0("nuc-1",tfix,rg);

EX.FX(x,tfix,rgw,sw)       =  ex0(x,tfix,rgw);
PRSV.FX(x,tfix,rgw,sw)      =  prsv0(x,tfix,rgw);
URSC.FX(x,tfix,rgw,sw)      =  ursc0(x,tfix,rgw);
RA.FX(x,tfix,rgw,sw)       =  ra0(x,tfix,rgw);

CG.FX("2000",sw) = 0;
CG.UP(pp,sw) = prsv0_w("coal") + ursc0_w("coal");
SC.FX(tfix,rg,sw) = 0;
SC.UP(pp,rg,sw) = sclim(rg);

*       Eliminate resource model (and trade) for aggregation levels not in use

EX.FX(x,pp,rgw,sw)$(not xmap(x,rgw)) = 0;
PRSV.FX(x,pp,rgw,sw)$(not xmap(x,rgw)) = 0;
URSC.FX(x,pp,rgw,sw)$(not xmap(x,rgw)) = 0;
RA.FX(x,pp,rgw,sw)$(not xmap(x,rgw)) = 0;
loop(x, NTX.FX(trd,pp,rg,sw)$(xtrd(x,trd) and (not xmap(x,rg))) = 0;);

*       Assign passenger transportation levels in the base year for those
*       regions in which this submodel is active:

PT.FX(pvt,tfix,trans(rg),sw) = pt0(pvt,tfix,rg);

* No new technologies can enter in 2010
PT.FX(pvt,"2010",trans(rg),sw)$(pt0(pvt,"2000",rg)=0) = 0;

* Define scenarios with no electric vehicles
$if set noev PT.FX("phev",pp,trans(rg),sw) = 0;
$if set noev PT.FX("elcv",pp,trans(rg),sw) = 0;

*       Fix passenger transportation levels in all other regions:

PT.FX(pvt,tp,ntrns(rg),sw) = 0;
PN.FX("evhcl",pp,rg,sw)$trans(rg) = 0;

*       Exogenous demand for transport:

T.FX(rg,pp,sw)$(trans(rg) and (not tpclay(rg))) = tref(pp,rg);
TC.FX(rg,pp,sw)$ntrns(rg) = 0;
T.M(rg,pp,sw) = 0;

* Net exports

NTX.FX("crude",tfix,rg,sw)  =  ntx0("crude",tfix,rg);
NTX.UP("crude",pp,rg,sw)$xmap("oil",rg)     =  oilx(pp,rg);
NTX.LO("crude",pp,rg,sw)$xmap("oil",rg)     =  - oilm(pp,rg);

NTX.FX("ngpl",tfix,rg,sw)  =  ntx0("ngpl",tfix,rg);
NTX.UP("ngpl",pp,rg,sw)     =  sum(rrg, ngplcap(pp,rrg,rg));
NTX.LO("ngpl",pp,rg,sw)     =  -sum(rrg, ngplcap(pp,rg,rrg));

NTX.FX("lng",tfix,rg,sw)  =  ntx0("lng",tfix,rg);
NTX.UP("lng",pp,rg,sw)     =  gasx(pp,rg);
NTX.LO("lng",pp,rg,sw)     =  -gasm(pp,rg);

EXPRT.FX(trd,tfix,rg,sw) = max(0, ntx0(trd,tfix,rg));
IMPRT.FX(trd,tfix,rg,sw) = -min(0, ntx0(trd,tfix,rg));

CLEV.FX(tfix,rg,sw) = clev0(tfix,rg);
CCEM.FX(tfix,rg,sw) = ccem0(tfix,rg);
CEMCCS.UP(tp,rg,sw) = cement(tp,rg);
CEMCCS.FX(tfix,rg,sw) = 0;
MLEV.FX(tfix,rg,sw) = mlev0(tfix,rg);
EM.FX("co2",tfix,sw) = carbem0(tfix);
EM.FX("ch4",tfix,sw) = sum(rg, mlev0(tfix,rg));

EM.FX("n2o",ctp,sw)        =  0;
EM.FX("slf",ctp,sw)        =  0;
EM.FX("llf",ctp,sw)        =  0;

ABATE.UP(ogg,abx,tp,rg,sw)  =   ablim(ogg,abx,tp,rg);
ABATE.FX("co2",abx,tp,rg,sw)   =  0;
ABATE.FX(ogg,abx,tfix,rg,sw) = abate0(ogg,abx,tfix,rg);

OTHCABT.UP(tp,sw) = bunker(tp) - bunker_min(tp);
OTHCABT.FX(tfix,sw) = 0;
$if set qty OTHCABT.FX(pp,sw) = othcabt_qty(pp,sw);

$if not set climate_fb NETLAND.FX(tp,sw) = 0;

*AFF.FX(tfix,ps,sw) = 0;
*CRLX.FX(tp,rg,sw) =  0;
*ORLX.FX(ogg,tp,rg,sw) =  0;

parameter        ceqxl(tp,rg)    Export limits on C-eq permit trading
                 ceqml(tp,rg)    Import limits on C-eq permit trading
                 ghgxl(trd,tp,rg)  Export limits on permit trading by gas
                 ghgml(trd,tp,rg)  Import limits on permit trading by gas;

ceqxl(pp,rg) = +inf;
ceqml(pp,rg) = +inf;

* No international trade in EMF 22 US scenario (except within Kyoto bubble)

$if set emf22us ceqxl(tp,rg)$(not kyt(rg)) = 0;  ceqml(tp,rg)$(not kyt(rg)) = 0;

loop(ghg, loop(xghg(ghg,trd), ghgxl(trd,pp,rg)$(not qtyout(pp,rg)) = +inf;
                              ghgml(trd,pp,rg)$(not qtyout(pp,rg)) = +inf;););

* Turn trade off after progressive target adoption
$if set fntrd ghgxl(trd,pp,nanb(rg)) = 0; ghgml(trd,pp,nanb(rg)) = 0;

$if %tenfx%==yes $goto skipkyoto
$if %kyoto%==no $goto skipkyoto

*       If 2010 is a decision year, a constraint can be
*       applied to represent Kyoto targets for [EUR, JAPAN, CANZ]:

set      kyt(rg);
kyt(rg)$(oecd(rg) and not sameas(rg,"usa")) = yes;

*       Abatement in all gases may be applied toward Kyoto target (using GWPs);
*       Trade is allowed only inside the Kyoto bubble (i.e. no CDM)

equation kytprot(rg,sw)  Kyoto Protocol for 2010 emissions;

kytprot(kyt(rg),sw)..
         CLEV("2010",rg,sw) - 0.001 * sum((ogg,abx), ABATE(ogg,abx,"2010",rg,sw)) =l=
                 kyoto(rg) - NTX("crt-ceq","2010",rg,sw);

ceqxl("2010",noecd) = 0;
ceqml("2010",noecd) = 0;
ceqxl("2010","usa") = 0;
ceqml("2010","usa") = 0;

$label skipkyoto

*ceqxl(tp,rg) = 0;
*ceqml(tp,rg) = 0;

loop(ghg, loop(xghg(ghg,trd), NTX.FX(trd,tfix,rg,sw) = 0;
                              NTX.UP(trd,pp,rg,sw)   =   ghgxl(trd,pp,rg);
                              NTX.LO(trd,pp,rg,sw)   = - ghgml(trd,pp,rg);););
NTX.FX("crt-ceq",tfix,rg,sw)     = 0;
NTX.UP("crt-ceq",pp,rg,sw)     =    ceqxl(pp,rg);
NTX.LO("crt-ceq",pp,rg,sw)     =  - ceqml(pp,rg);

NBC.FX(tp,rg,sw)$(not tnbc(tp,rg)) = 0;
$if set nobrw CBC.LO(tp,"usa",sw) = 0;

$if %mktdam%==no  MD.FX(rg,ctp,sw) = 0;
$if %nmktdam%==no ELF.FX(rg,ctp,sw) = 1;
$if %nmktdam%==no  udf(pp,rg) = udf(pp,rg)/sum(tp$pp(tp),udf(tp,rg));
$if not set mu RFMAX.FX = 0;
$if not set tau TPMAX.FX = 0;
$if not set zeta CONMAX.FX = 0;

$ontext
$if set rfktrg $setglobal rfktrg %rfktrg%
$if set tmptrg $setglobal tmptrg %tmptrg%
$if set mu $setglobal mu %mu%
$if set tau $setglobal tau %tau%
$if set zeta $setglobal zeta %zeta%
$offtext

* When pfix has been turned on, energy variables are held at reference case for certain regions
$if exist pfix_%pfix%.gms $include pfix.gms
$if %pfix%==yes $include pfix.gms

* * * * * * * * * * * * * * Model Execution * * * * * * * * * * * * *

MODEL M7 /ALL/;

OPTION NLP = CONOPT3;
M7.OPTFILE  = 1;
OPTION SOLPRINT  =    on;
OPTION LIMROW    =      0;
OPTION LIMCOL    =      0;
OPTION ITERLIM   =  50000;
OPTION RESLIM    =  50000;

*        Initialize iterative parameters
$if not defined kterm loop(tlast, kterm(rg,sw) = K.L(rg,tlast,sw)*(1+grow(tlast,rg))**10;);
$if set nash ghgtax(ghg,pp,rg,sw)$(ghgtax(ghg,pp,rg,sw) eq 0) = -nashsbd(ghg,pp,rg,sw);
ghgtax(ghg,psub(tp,rg),sw) = -ghgsbd(ghg,tp,rg,sw);
$if set gradjoin ghgtax(ghg,pp(tp),nanb(rg),sw)$(not psub(tp,rg) and not pfix(tp,rg)) = -ghgsbd(ghg,tp,rg,sw);


$if set climate parameter emitlog Iteration log of emission paths;
$if set climate parameter lambda  Post-terminal adjustment step length (full step initially) /1/;
$if %lifeiter%==yes parameter lifeiter Iteration log of CH4 and N2O lifetimes;
$if %tmfbiter%==yes parameter tmfbiter Iteration log of temperature feedback;

* Initialize Nash coalition set
$if set nash $include marginal_e.gms

*       Generate a status report in the title bar if we are operating
*       on Windows NT/2K/XP:

$if %system.filesys%==MSNT file title_ /'title.cmd'/; title_.nd=0; title_.nw=0;
$set updatetitle "putclose title_ '@title Solving %basisout%, iter ',ord(iter),' of ',card(iter),'. Start time: %system.time%, Current time: %time% -- Ctrl-S to pause'/; execute 'title.cmd';"


$ontext
* * * * * Check climate variables * * * * * *
execute_unload 'climtest_in.gdx', TOTEM.L, CO2.L, S.L, SMINUS.L;
$include rpt_climate.gms
execute_unload 'climtest_out.gdx', TOTEM.L, CO2.L, S.L, SMINUS.L;
$offtext


parameter       encompare       Comparison of E and N with reference path;
encompare("ref","e",tp,rg)  = eref(tp,rg);
encompare("ref","n",tp,rg)  = nref(tp,rg);
encompare("ref","en",tp,rg) = enref(tp,rg);
encompare("ref","nn",tp,rg) = nnref(tp,rg);

*.M7.holdfixed = yes;
*.OPTION NLP=EXAMINER;
*.SOLVE M7 MAXIMIZING NWEL USING NLP;
*.$exit

loop(iter,

$if defined title_ %updatetitle%

*$if %nmktdam%==yes U.FX(rg,ctp,sw)$((not tp(ctp)) and tpclay(rg)) = (U.L(rg,"2100",sw)+(C.L(rg,"2100",sw)/(C.L(rg,"2100",sw) + I.L(rg,"2100",sw))) * MD.L(rg,"2100",sw)) * 1.01**(yr(ctp) - yr("2100"))-MD.L(rg,ctp,sw);
*$if %nmktdam%==yes C.FX(rg,ctp,sw)$((not tp(ctp)) and (not tpclay(rg))) = (C.L(rg,"2100",sw)*(1 + MD.L(rg,"2100",sw)/(C.L(rg,"2100",sw) + I.L(rg,"2100",sw)))) * 1.01**(yr(ctp) - yr("2100"))-MD.L(rg,ctp,sw);

        nwtitr(iter,rg)       = nwt(rg);

*.        M7.holdfixed = yes;
        SOLVE M7 MAXIMIZING NWEL USING NLP;
        abort$(M7.solvestat=5) "Model is infeasible.";

$if not set climate $goto iterref

*        Calculate subsidies for all GHGs for non/partial participants

        ghgsbd(ghg,psub(tp,rg),sw) = -1000 * totemit.m(ghg,tp,sw)/abs(trdbal.M("nmr",tp,sw));
$if %gradjoin%==ppc ghgsbd(ghg,pp(tp),nanb(rg),sw)$(not psub(tp,rg) and not pfix(tp,rg)) = (1 - ppc_rel(tp,rg)) * (-1000 * totemit.m(ghg,tp,sw)/abs(trdbal.M("nmr",tp,sw)));
$if %gradjoin%==emf ghgsbd(ghg,pp(tp),nanb(rg),sw)$(not psub(tp,rg) and not pfix(tp,rg)) = grdjnemf(tp,rg) * (-1000 * ((totemit.m(ghg,tp,sw)/abs(trdbal.M("nmr",tp,sw))) - (totemit.m(ghg,"2020",sw)/abs(trdbal.M("nmr","2020",sw)))));


$if set nash $include marginal_e.gms

*       Project post-terminal emissions (geometric decline):

        emitlog(ghg,ctp,iter,sw) = TOTEM.L(ghg,ctp,sw);
        loop(ctp$(not tp(ctp)), TOTEM.FX(ghg,ctp,sw) = (1-lambda) * TOTEM.L(ghg,ctp,sw) + lambda * TOTEM.L(ghg,ctp-1,sw) * 0.9;);
* min(1,TOTEM.L(ghg,ctp-1,sw)/TOTEM.L(ghg,ctp-2,sw)));
*       When stabilizing temperature, need to adjust first attempt at post-terminal emissions
$if set tmptrg if(ord(iter) eq 1, TOTEM.FX(ghg,ctp,sw)$(not tp(ctp)) = 0.25 * TOTEM.L(ghg,ctp,sw););

* Include extended climate model and set step length
         cpp(ctp)$(not ctfix(ctp)) = yes;
         ecp(cpp) = yes; loop(tfix, ecp(ctp)$sameas(ctp,tfix) = no;);
         lambda = 0.75;

* If time-variants are updated manually, include iteration step
$if %lifeiter%==yes $include lifeiter.gms
$if %tmfbiter%==yes $include tmfbiter.gms

$label iterref

        loop((p0,ssw)$(ord(ssw)=1),
          pvpi(trd,pp,sw) = sum(swm(sw,sow,pp),
                abs(TRDBAL.M(trd,pp,sow)/TRDBAL.M("nmr",p0,ssw)));
          pvpt(pp,trans(rg),sw) =  sum(swm(sw,sow,pp),
                abs(SUPTRANS.M(rg,pp,sow)/TRDBAL.M("nmr",p0,ssw)));
        );

        nwt(rg) = sum((pp,sw), sum(swm(sw,sow,pp),
                (pvpi("nmr",pp,sow)*C.L(rg,pp,sow) + (pvpt(pp,rg,sw)*T.L(rg,pp,sw))$trans(rg))
*  ??           *ELF.L(rg,pp,sow)
                + sum(trd, pvpi(trd,pp,sow)*NTX.L(trd,pp,rg,sow))))
$if set mu + sh("2020",rg) * mu * RFMAX.L
$if set tau + sh("2020",rg) * tau * TPMAX.L
$if set zeta + sh("2020",rg) * zeta * CONMAX.L
        ;
        nwt(rg) = nwt(rg) / sum(rrg, nwt(rrg));

*       Update xshr for one-world resources

        xshr(x,tp,rg,sw) = 0;
        xshr(x,tp,rg,sw)$sum(rrg, TPE.L(x,tp,rrg,sw)) =
                TPE.L(x,tp,rg,sw) / sum(rrg, TPE.L(x,tp,rrg,sw));

*       Target the terminal capital stock:

        loop(tlast, kterm(rg,sw) = K.L(rg,tlast,sw)*(1+grow(tlast,rg))**10; );

*       Calibrate a disutility of nuclear power such that when
*       nuclear output represents benchmark share of electric supply, the
*       marginal risk premium is wtp0 (10 mills per kwh in OECD):

        rsktax(rsk,pp,rg,sw) = (wtp0(rsk,rg)/2) * (1/ (refshr(rsk,rg)*E.L(rg,pp,sw)));

*       Include adjustments for emissions subsidies representing partial participation
$if set nash ghgtax(ghg,pp,rg,sw) = -nashsbd(ghg,pp,rg,sw);
        ghgtax(ghg,psub(tp,rg),sw) = -ghgsbd(ghg,tp,rg,sw);
$if set gradjoin ghgtax("co2",pp,nanb(rg),sw)$(not psub(pp,rg) and not pfix(pp,rg)) = -ghgsbd("co2",pp,rg,sw);

        ghgitr(ghg,tp,rg,iter,sw) = ghgtax(ghg,tp,rg,sw);

        taxrev(pp(tp),rg,sw) = sum(swm(sw,sow,tp),
                  xntax(rg)*N.L(rg,tp,sow)
                + xetax(rg)*E.L(rg,tp,sow)
                + ghgtax("co2",tp,rg,sow)*(CLEV.L(tp,rg,sow) + CCEM.L(tp,rg,sow))
*               + ghgtax("co2",tp,rg,sow)*sh(tp,rg)*AFFEMIT.L(tp,sw)
                + ghgtax("ch4",tp,rg,sow) * MLEV.L(tp,rg,sow)
                + sum(ogg(ghg), ghgtax(ghg,tp,rg,sow) * (OTHEMIT.L(ghg,tp,rg,sw) - sum(ogs, oggmin(ghg,ogs,tp,rg))))
                + rsktax("n",tp,rg,sow)*sqr(sum(nuc, PE.L(nuc,tp,rg,sow)))
                + rsktax("c",tp,rg,sow)*sqr(sum(ccs, PE.L(ccs,tp,rg,sow)))
                + rsktax("p",tp,rg,sow)*sqr(PE.L("nuc-adv",tp,rg,sow)));

        taxitr(tp,rg,iter,sw) = taxrev(tp,rg,sw);

        iprev(tp,rg,sw) = sum(trd, ipprem(trd,tp,rg)* EXPRT.L(trd,tp,rg,sw));

$if not set recal $goto enditer
* * * * * Recalibrate to exogenous baseline energy demand path * * * * * * * *

        loop(sw,
*        Observe prices for energy and levels for capital and output
          peref(pp,rg) = (newelec.m(rg,pp,sw) / newprod.m(rg,pp,sw))$%gross% +
                         (newelec.m(rg,pp,sw) / newprod.m(rg,pp,sw))$%newvint%;
          pnref(pp,rg) = (newnon.m(rg,pp,sw) / newprod.m(rg,pp,sw))$%gross% +
                         (newnon.m(rg,pp,sw) / newprod.m(rg,pp,sw))$%newvint%;
          yref(pp,rg) = Y.L(rg,pp,sw);
          ynref(pp,rg) = YN.L(rg,pp,sw);
          kref(pp,rg) = K.L(rg,pp,sw);
          knref(pp,rg) = KN.L(rg,pp,sw);
*       Recalculate value shares based on observations and target energy path
          theta(pp,rg) =
                 ((peref(pp,rg)*eref(pp,rg)+pnref(pp,rg)*nref(pp,rg)) / yref(pp,rg))$%gross% +
                 ((peref(pp,rg)*enref(pp,rg)+pnref(pp,rg)*nnref(pp,rg)) / ynref(pp,rg))$%newvint%;
          elvs(pp,rg) =
                 ((peref(pp,rg)*eref(pp,rg)) / (peref(pp,rg)*eref(pp,rg) + pnref(pp,rg)*nref(pp,rg)))$%gross% +
                 ((peref(pp,rg)*enref(pp,rg)) / (peref(pp,rg)*enref(pp,rg) + pnref(pp,rg)*nnref(pp,rg)))$%newvint%;
        );

$label enditer

        encompare(iter,"e",tp,rg) = sum(sw,E.L(rg,tp,sw))/card(sw);
        encompare(iter,"n",tp,rg) = sum(sw,N.L(rg,tp,sw))/card(sw);
        encompare(iter,"en",tp,rg) = sum(sw,EN.L(rg,tp,sw))/card(sw);
        encompare(iter,"nn",tp,rg) = sum(sw,NN.L(rg,tp,sw))/card(sw);

);

$if %iter%==it0 $exit

$if set climate execute_unload 'emitlog',emitlog;
option decimals = 6;
DISPLAY nwtitr;

*       We have completed the baseline solution.  Do one more solve
*       to provide equation listing, column listing, solution listing,
*       shadow prices and a basis file:

OPTION SOLPRINT  =      on;
OPTION LIMROW    =     100;
OPTION LIMCOL    =     100;
*.M7.holdfixed = no;
M7.savepoint = 1;
SOLVE M7 MAXIMIZING NWEL USING NLP;

* Show main outputs
parameter primal(tp,rg,*);
primal(pp,rg,"rlzgdp") = C.L(rg,pp,"ref") +  I.L(rg,pp,"ref") + TC.L(rg,pp,"ref")
        + sum(trd, NTX.L(trd,pp,rg,"ref") * abs(trdbal.M(trd,pp,"ref")/trdbal.M("nmr",pp,"ref"))) ;
primal(tfix,rg,"rlzgdp") = gdp(tfix,rg);
primal(tfix,rg,"rlzppc") = ppc(tfix,rg);
loop(pp(tp), primal(tp,rg,"rlzppc") = primal(tp-1,rg,"rlzppc") * (primal(tp,rg,"rlzgdp")/primal(tp-1,rg,"rlzgdp")) /
                                                                 (pop(tp,rg) / pop(tp-1,rg)););
primal(tp,rg,"tpe") = sum(v, sum(etv, htrt_v(etv,v,tp,rg) * PEV.L(etv,v,tp,rg,"ref")))
                    + sum(et$(not etv(et)), 10 * PE.L(et,tp,rg,"ref"))
                    + N.L(rg,tp,"ref") + PR.L("synf",tp,rg,"ref")*syntpe
                    + sum(pvt, sum(ntf(tf), epvkt(tf,pvt,tp,rg) * PT.L(pvt,tp,rg,"ref")));
primal(tp,rg,"pctpe") = primal(tp,rg,"tpe") / pop(tp,rg);
primal(tp,rg,"co2") = CLEV.L(tp,rg,"ref") + CCEM.L(tp,rg,"ref");

parameter htrt(et,tp,rg) average heat rate;
htrt(etf,tp,rg)$PE.L(etf,tp,rg,"ref") = sum(vtp(v,tp),
                 htrt_v(etf,v,tp,rg) * PEV.L(etf,v,tp,rg,"ref")) / PE.L(etf,tp,rg,"ref");

* Show prices
parameter dual(tp,rg,*);
dual(pp,rg,"elec")   = 1000*supelec.M(rg,pp,"ref")/cc.M(rg,pp,"ref");
dual(pp,rg,"nele")   = 1000*supnon.M(rg,pp,"ref")/cc.M(rg,pp,"ref");
dual(pp,rg,"liq")   = 1000*supliq.M(rg,pp,"ref")/cc.M(rg,pp,"ref");
$if not set woil dual(pp,rg,"oil")   = 1000*6.0*supx.M("oil",rg,pp,"ref")/cc.M(rg,pp,"ref");
$if set woil dual(pp,rg,"oil")   = 1000*6.0*wsupx.M("oil",pp,"ref")/trdbal.M("nmr",pp,"ref");
$if not set wgas dual(pp,rg,"gas")   = 1000*supx.M("gas",rg,pp,"ref")/cc.M(rg,pp,"ref");
$if set wgas dual(pp,rg,"gas")   = 1000*wsupx.M("gas",pp,"ref")/trdbal.M("nmr",pp,"ref");
dual(pp,rg,"uranium") = 1000* (500/2.2) * wsupx.M("ur",pp,"ref")/trdbal.M("nmr",pp,"ref");
dual(pp,rg,"nuc-elec") = ecst_fuel("nuc-1",pp,rg) + 1000 * wsupx.M("ur",pp,"ref")/trdbal.M("nmr",pp,"ref");
dual(pp,rg,"trans")   = 1000*suptrans.M(rg,pp,"ref")/cc.M(rg,pp,"ref");

dual(pp,rg,"carbon")$(not pfix(pp,rg)) = ghgtax("co2",pp,rg,"ref") +
        1000 * (-wcardf.m(pp,"ref") + clevbd.m(rg,pp,"ref") + ceqbd.m(rg,pp,"ref"))
        / abs(trdbal.m("nmr",pp,"ref"));

dual(pp,rg,"ceq") = 1000 * (ceqbd.m(rg,pp,"ref"))/abs(trdbal.m("nmr",pp,"ref"));

parameter premium(rsk,tp,rg)     Risk premium applied at the margin;
premium("n",tp,rg) = 2 * rsktax("n",tp,rg,"ref") * sum(nuc, PE.L(nuc,tp,rg,"ref"));
premium("c",tp,rg) = 2 * rsktax("c",tp,rg,"ref") * sum(ccs, PE.L(ccs,tp,rg,"ref"));
premium("p",tp,rg) = 2 * rsktax("p",tp,rg,"ref") * PE.L("nuc-adv",tp,rg,"ref");

parameter transrpt(tp,rg,*)      Transportation Energy report;
transrpt(tp,rg,"TKWh-final") = ( eveff(tp,rg) * PN.L("evhcl",tp,rg,"ref") )$ntrns(rg)
        + (1/3.6) * sum(pvt, epvkt("elec",pvt,tp,rg) * PT.L(pvt,tp,rg,"ref"))$trans(rg);
transrpt(tp,rg,"E-EJ-final") = 3.6 * ( eveff(tp,rg) * PN.L("evhcl",tp,rg,"ref") )$ntrns(rg)
                  + sum(pvt, epvkt("elec",pvt,tp,rg) * PT.L(pvt,tp,rg,"ref"))$trans(rg);
transrpt(tp,rg,"N-EJ-final") = sum(ntf, sum(pvt, epvkt(ntf,pvt,tp,rg) * PT.L(pvt,tp,rg,"ref")));
transrpt(tp,rg,"N-EJ-prim") = sum(ntf, sum(pvt, epvkt(ntf,pvt,tp,rg) * PT.L(pvt,tp,rg,"ref"))) +
                              0.66 * sum(pvt, epvkt("lqtr",pvt,tp,rg) * PT.L(pvt,tp,rg,"ref")) *
                              PR.L("synf",tp,rg,"ref") / sum(lf, PR.L(lf,tp,rg,"ref"));

* Define and assign climate variables for reference runs
$if not set climate $include rpt_climate.gms
$if not set climate $goto nonash
$if %mktdam%==no  MD.L(rg,ctp,sw)  = mdmfac(rg) * gdp(ctp,rg) * ATP.L(ctp,sw)/damage("reftemp",sw);
$if %nmktdam%==no ELF.L(rg,ctp,sw) = (1 - damage("refwtp",sw) * (ATP.L(ctp,sw)/damage("reftemp",sw))**2)**hsx(rg,ctp,sw);
$label nonash
$if not set nash parameter mrgben_rg(ghg,tp,rg); mrgben_rg(ghg,tp,rg) = 0;

execute_unload 'm7_nwt',nwt, taxrev, rsktax, ghgsbd, nashsbd, iprev, kterm;
$if not set study $set study
$if set basisout execute 'copy m7_p.gdx %basdir%%basisout%.gdx';
$if set basisout execute 'copy m7_nwt.gdx %basdir%%basisout%_nwt.gdx';

execute 'if exist title.cmd del title.cmd';
execute 'copy m7_p.gdx %basdir%m7_p.gdx';
execute 'copy m7_nwt.gdx %basdir%m7_nwt.gdx';
execute '=gdx2xls m7_p.gdx m7_p.xlsx';
execute '=gdx2xls m7_nwt.gdx m7_nwt.xlsx';