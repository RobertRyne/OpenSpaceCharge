c***********************************************************************
c
c Bessel functions of the first kind.
c D.T.Abell, Tech-X Corp., November 2003.
c
c Updated December 2005---improved rational approximations, obtained by
c comparison with Mathematica's high-precision Bessel function values.
c
c***********************************************************************
      subroutine dbessj0(x,bj0)
c Bessel function of the first kind, order zero: bj0 = J_0(x).
      implicit none
      double precision x,bj0
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      double precision Pi,Pid4
      double precision ax,t,pp,qq
c
c coefficients for rational approximations
      double precision rabj0an(11),rabj0ad(11),rabj0bn(11),rabj0bd(11)
      double precision rabj0cn(11),rabj0cd(11),rabj0dn(11),rabj0dd(11)
      double precision rabj0en(11),rabj0ed(11),rabj0fn( 7),rabj0fd( 7)
      data rabj0an /
     &+1.0000000000000000000050179699999091D+00,
     &-8.8230193349303631019069286292273592D-02,
     &-2.2778858233747603162266489897356056D-01,
     &+2.0341718595179019764036039758709466D-02,
     &+1.0311930796055483981361241333992624D-02,
     &-9.6558198466558060525000633951901562D-04,
     &-1.4534545971264389949915743465588711D-04,
     &+1.5385739184044285945395365895073623D-05,
     &+5.0024280550288686235912958144335425D-07,
     &-8.1837015732913880463394834114513691D-08,
     &+1.8402156332042362862773578704341447D-09/
      data rabj0ad /
     &+1.0000000000000000000000000000000000D+00,
     &-8.8230193349303629919583311940311153D-02,
     &+2.2211417662523928303822343473518375D-02,
     &-1.7158297421463074264424446502031107D-03,
     &+2.3978521168201685670512690979806517D-04,
     &-1.5942649098416680654142008778904110D-05,
     &+1.5752199434140609397499317066559390D-06,
     &-8.4437977811612729526150890484969979D-08,
     &+6.0918476815384520168977485010843715D-09,
     &-2.1081214567148162491668339094873657D-10,
     &+9.5689922998775932286006136295490244D-12/
      data rabj0bn /
     &+0.9999947618451229408378453732245815D+00,
     &-2.1079067749100990034907737987252581D-01,
     &-2.1153610058761303594734319673660807D-01,
     &+4.7914960041618770111290430390100304D-02,
     &+6.4925682382712033512704530646438245D-03,
     &-2.1181715676820241543828289838391954D-03,
     &+3.5262828624115159195904623163278898D-05,
     &+2.8964640803552355496585442664205909D-05,
     &-3.1467459340839512757355221349774376D-06,
     &+1.3056612961826745813193969076704308D-07,
     &-1.9723613039425622244289372840306964D-09/
      data rabj0bd /
     &+1.0000000000000000000000000000000000D+00,
     &-2.1080954250115548577070508093772022D-01,
     &+3.8495956449617845060277310851005345D-02,
     &-4.8215559602633320002172981719561965D-03,
     &+5.1704518623798707902299649434725156D-04,
     &-4.3816500419579022209215618337287690D-05,
     &+3.0961318132620954195784080539882816D-06,
     &-1.6985374380267480945014770381800044D-07,
     &+7.1683418472403236347208318913158728D-09,
     &-2.0143582141730917107139467960886491D-10,
     &+3.2401851766762590203597072719481002D-12/
      data rabj0cn /
     &+4.8896470117798136933960119243767237D-02,
     &-1.9199232842235873850450795814434952D-03,
     &-5.7862889184410290524315292285640746D-02,
     &+3.9102508961576130540892154545854854D-02,
     &-1.1494590155653841967222931097121521D-02,
     &+1.8803934283457853100848553874229646D-03,
     &-1.8640431669347982491495399716692452D-04,
     &+1.1474374355794768767774947971846314D-05,
     &-4.2933509072418537179570698396689563D-07,
     &+8.9528106522417889070362524823241969D-09,
     &-7.9882819654637042538953026938994769D-11/
      data rabj0cd /
     &+1.0000000000000000000000000000000000D+00,
     &-1.1329555327217814998913315706968691D+00,
     &+3.2591109585817699371371552353553638D-01,
     &-5.7426238664400627746092218763981592D-02,
     &+6.9805158757947823161383120930805137D-03,
     &-6.1938100100437024756343054169434183D-04,
     &+4.0427188509225439620395057691850240D-05,
     &-1.9137212931564333474744008597798173D-06,
     &+6.2684456139853366202327362897514514D-08,
     &-1.2847022551902166826121217722839446D-09,
     &+1.2563046316660055826508853359682238D-11/
      data rabj0dn /
     &-1.2855106884992982411556047079499048D-02,
     &+7.3122100400657976523877937106051494D-02,
     &-5.3227204991088052536112979278795060D-02,
     &+1.7018649354405257105380788618114423D-02,
     &-3.0419674749227155566039454550406450D-03,
     &+3.3475526134474814407534241200383970D-04,
     &-2.3600739767436000858506331588751560D-05,
     &+1.0701657071023291919125855988389812D-06,
     &-3.0207072100726981004102397058083471D-08,
     &+4.8321242798880922243680861245207648D-10,
     &-3.3474933906957231905826150510301703D-12/
      data rabj0dd /
     &+1.0000000000000000000000000000000000D+00,
     &-3.2944072183253369366766487533934724D-01,
     &+6.0965775940009661518427458794047474D-02,
     &-7.5921317989016470654981426275636095D-03,
     &+6.8802934300040013792164096446792849D-04,
     &-4.6422034672029920649858016846729019D-05,
     &+2.3380523557789858162629933228290598D-06,
     &-8.6093073725985170178907421779610271D-08,
     &+2.2141159560680130851824587967869377D-09,
     &-3.5909539384246578384084604881458476D-11,
     &+2.8364122289958662381014742871402263D-13/
      data rabj0en /
     &+3.5224646241387634580428630255239534D-01,
     &-2.7856689678561711255883463632505246D-01,
     &+9.3729467733125313886079046699370299D-02,
     &-1.7776325167929891943907611157376352D-02,
     &+2.1156292810051315475968635296554389D-03,
     &-1.6585801993726409351607592862093892D-04,
     &+8.7086646940913853360940919270539488D-06,
     &-3.0345030581942515926188590998173107D-07,
     &+6.7356307310329723356037410610211086D-09,
     &-8.6230544824987840642736370362678861D-11,
     &+4.8463254940494631096363486040475736D-13/
      data rabj0ed /
     &+1.0000000000000000000000000000000000D+00,
     &-3.9437749159126193952363812178419292D-01,
     &+7.7581865230995871051567171336090776D-02,
     &-9.6726853312945266401280401225401609D-03,
     &+8.3352527419529980773115111930086124D-04,
     &-5.1371247731535664030011300083139149D-05,
     &+2.2801094464620980467671861325928996D-06,
     &-7.1714183814137341630980676930651847D-08,
     &+1.5275616049573480145724295285266087D-09,
     &-1.9906555330977615027102660888016815D-11,
     &+1.2106511741132039245923717952671393D-13/
      data rabj0fn /
     &-9.4845363290573921489820379821276979D-02,
     &+2.8600548026603177758882128713536768D-02,
     &-3.5476668104961868720495540952881358D-03,
     &+2.3179029870614323398028994082351304D-04,
     &-8.4165980090902303831222302065086657D-06,
     &+1.6111527211980360198505293917220034D-07,
     &-1.2708319802262517343077780038561572D-09/
      data rabj0fd /
     &+1.0000000000000000000000000000000000D+00,
     &-2.6409647927463548120030732742120998D-01,
     &+2.9860012401515139424647801546743694D-02,
     &-1.8401803786971692142522329594649180D-03,
     &+6.5061344201723819188622837380416900D-05,
     &-1.2504130549754413803085326308427806D-06,
     &+1.0221043546656880551547140513637591D-08/
c
c compute constants
      Pi=2.d0*dasin(1.d0)
      Pid4=0.25d0*Pi
c
c J_0(x) is even in x, so we use |x|
      ax=dabs(x)
c
c use approximation appropriate to value of x:
c   range (0,4): rational approximation of J_0(x)
c   range [4,8): rational approximation of J_0(x)
c   range [8,12): rational approximation of
c                    J_0(x)-sqrt(2/(pi*x))*cos(x-pi/4)
c   range [12,15): rational approximation of
c                    J_0(x)-sqrt(2/(pi*x))*cos(x-pi/4)
c   range [15,20): rational approximation of
c                    J_0(x)-sqrt(2/(pi*x))*cos(x-pi/4)
c   range [20,21): rational approximation of
c                    J_0(x)-sqrt(2/(pi*x))*cos(x-pi/4)
c   range [21,infinity): Hankel's asymptotic expansion
c
      if (ax.eq.0.d0) then
        bj0=1.d0
      else if (ax.lt.4.d0) then
        call ratappr(ax,rabj0an,11,rabj0ad,11,bj0)
      else if (ax.lt.8.d0) then
        call ratappr(ax,rabj0bn,11,rabj0bd,11,bj0)
      else if (ax.lt.12.d0) then
        call ratappr(ax,rabj0cn,11,rabj0cd,11,bj0)
        bj0=bj0+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pid4)
      else if (ax.lt.15.d0) then
        call ratappr(ax,rabj0dn,11,rabj0dd,11,bj0)
        bj0=bj0+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pid4)
      else if (ax.lt.20.d0) then
        call ratappr(ax,rabj0en,11,rabj0ed,11,bj0)
        bj0=bj0+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pid4)
      else if (ax.lt.21.d0) then
        call ratappr(ax,rabj0fn,7,rabj0fd,7,bj0)
        bj0=bj0+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pid4)
      else
        t=ax-Pid4
        call HankelPQ(0,ax,pp,qq)
        bj0=dsqrt(2.d0/(Pi*ax))*(pp*dcos(t)-qq*dsin(t))
      end if
c
      return
      end
c
c***********************************************************************
      subroutine dbessj1(x,bj1)
c Bessel function of the first kind, order one: bj1 = J_1(x).
      implicit none
      double precision x,bj1
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      double precision Pi,Pi3d4
      double precision ax,t,pp,qq,sgn
c
c coefficients for rational approximations
      double precision rabj1an(11),rabj1ad(11),rabj1bn(11),rabj1bd(11)
      double precision rabj1cn(11),rabj1cd(11),rabj1dn(11),rabj1dd(11)
      double precision rabj1en(11),rabj1ed(11),rabj1fn( 7),rabj1fd( 7)
      data rabj1an /
     &+1.3024604716812625136575985066550351D-21,
     &+4.9999999999999999971769784794883053D-01,
     &-2.7021475700860125161066293561710086D-02,
     &-5.3285436029298919304319635517409686D-02,
     &+2.9036486061909128531536930260765731D-03,
     &+1.5368020011789656151354698315581584D-03,
     &-8.5478332785495533785549719799392937D-05,
     &-1.6331599044080528432091645493034297D-05,
     &+9.4299064900371993501932762368482883D-07,
     &+5.9154546709636257637272914412136624D-08,
     &-3.6255546773750573374562132519289283D-09/
      data rabj1ad /
     &+1.0000000000000000000000000000000000D+00,
     &-5.4042951401720270639311664158490582D-02,
     &+1.8429127941402451229391111842933448D-02,
     &-9.4807171283539186854248011429003618D-04,
     &+1.6891166170994095677427891725156792D-04,
     &-7.9919244887708660985328880211035730D-06,
     &+9.7266277426087547516841459653563922D-07,
     &-3.9171386712324223695277101066676514D-08,
     &+3.4953779352137137722880097405620162D-09,
     &-9.4952410837245606026644100187155292D-11,
     &+5.8322537268273354279510519335860106D-12/
      data rabj1bn /
     &-9.6126582523756719172017705951900558D-06,
     &+5.0003360984690369391627958263960859D-01,
     &-1.0135247312540184243931882056358513D-01,
     &-4.2699412675762871627685186040163278D-02,
     &+1.0044781754253024712239985871722026D-02,
     &+4.5619419956893940475029673398794519D-04,
     &-2.4129439668031979459455941679644064D-04,
     &+1.6016153674406627617877625760637880D-05,
     &+2.1027192253057843266329778129856985D-07,
     &-5.0707491857311754051124557929869265D-08,
     &+1.2617640048087146823185552850801499D-09/
      data rabj1bd /
     &+1.0000000000000000000000000000000000D+00,
     &-2.0259445599160107275639117333790023D-01,
     &+3.9487894299400709115590578713486066D-02,
     &-5.1537400250629675578459531550226451D-03,
     &+5.9728087755675845907468596531647125D-04,
     &-5.4444672786583878044335341357153240D-05,
     &+4.2045344283680687094812530071394440D-06,
     &-2.5132379473265406943787755252580057D-07,
     &+1.1653755299489430413597047842196229D-08,
     &-3.5807396467498635498582625262060225D-10,
     &+6.3908259113402294201974014102826536D-12/
      data rabj1cn /
     &+6.0371438405985520492273239641094956D-01,
     &-5.9666583821895112063164472546734882D-01,
     &+2.6919847622163680583156066524301151D-01,
     &-8.7672496487695919979136196550569105D-02,
     &+2.1837199551259834491158951946045927D-02,
     &-3.7370808859808720836608377144630291D-03,
     &+4.1191847569201912836238278879683956D-04,
     &-2.8438941172117570987518632348259140D-05,
     &+1.1844617605399366234082982326828033D-06,
     &-2.7169771614418824284498359366424714D-08,
     &+2.6343527972027239886711798412892213D-10/
      data rabj1cd /
     &+1.0000000000000000000000000000000000D+00,
     &-1.6514728561739773326127601166782867D-01,
     &+6.5774300517547309014897346964600873D-02,
     &-1.5049813610294697551699714956373893D-02,
     &+2.4206604468666677360733935930230549D-03,
     &-2.6756346375625205474234901528281541D-04,
     &+2.1308038490017812103019561719321294D-05,
     &-1.2019599969650507783017925415454145D-06,
     &+4.6774324464657799887986331519504851D-08,
     &-1.1350580273767823974385743078813233D-09,
     &+1.3694521779905435685524410856631591D-11/
      data rabj1dn /
     &+4.2467609762315442400693029865795040D-01,
     &-4.3083142135486853345320553054798065D-01,
     &+1.6788517206101123528239344102761124D-01,
     &-3.3796177537444951633168287614832691D-02,
     &+3.8795047097080453256087993331777921D-03,
     &-2.5515818965006319689426544517178161D-04,
     &+8.3162429845444497107382145316229273D-06,
     &-3.1934959167383968458778078011436116D-09,
     &-9.3164847018358148213745655834193359D-09,
     &+2.8746610644386035220167446103613304D-10,
     &-2.8944287034616799451929693172577113D-12/
      data rabj1dd /
     &+1.0000000000000000000000000000000000D+00,
     &-3.6698450390768612683634161865311059D-01,
     &+7.4884143212216256039432139181885445D-02,
     &-1.0221916096061684403865100654999128D-02,
     &+1.0074562131863317846318355992966287D-03,
     &-7.3374157538561828451182655688029948D-05,
     &+3.9570224775017050301357112675936376D-06,
     &-1.5479721151459859100101554063875132D-07,
     &+4.1934519331039533720220137503174087D-09,
     &-7.1015545891050449440380478375745093D-11,
     &+5.7826067387808782775528298874601191D-13/
      data rabj1en /
     &-1.0727178100039883419848523490500134D+00,
     &+7.0928196569351176003941921501759760D-01,
     &-2.0244112729254919784586618692965949D-01,
     &+3.2919122179062670197381146685815259D-02,
     &-3.3846375280985925074645334995793181D-03,
     &+2.3035146427123168406922044131847732D-04,
     &-1.0526097486773898355270589440366091D-05,
     &+3.1928037549687904672248561941454213D-07,
     &-6.1569502122571540023308673784804519D-09,
     &+6.8176863727717420052494188810057293D-11,
     &-3.2903939595120148107005430367355315D-13/
      data rabj1ed /
     &+1.0000000000000000000000000000000000D+00,
     &-4.1294556742055131913125357418100254D-01,
     &+8.2841362688762866622842985663978347D-02,
     &-1.0382391872762524304804231048703553D-02,
     &+8.9012846103463523752825990646025319D-04,
     &-5.4163294972308871532368229256156985D-05,
     &+2.3586205292066095764629097580766725D-06,
     &-7.2385925159801043526475622512870661D-08,
     &+1.4965837356032315825582092468612243D-09,
     &-1.8825359032178275136083152292764955D-11,
     &+1.0969526078193747331991670408293388D-13/
      data rabj1fn /
     &-9.8075828607639007338780275479940338D-03,
     &+6.2042628485312133436957794665065060D-03,
     &-1.1711582212026499092489447079398142D-03,
     &+1.0277604893169514981220012453729818D-04,
     &-4.6894343161978861506416373649349031D-06,
     &+1.0823475910761513088113165275317871D-07,
     &-1.0009149629230385757860595107751339D-09/
      data rabj1fd /
     &+1.0000000000000000000000000000000000D+00,
     &-2.7630107315667408326367307184433078D-01,
     &+3.2444450997756793186658519895099244D-02,
     &-2.0654832225935412940479575463019885D-03,
     &+7.5094528449918628311451893661611736D-05,
     &-1.4778957595810851829028146021523979D-06,
     &+1.2314771882129871146918134466632580D-08/
c
c compute constants
      Pi=2.d0*dasin(1.d0)
      Pi3d4=0.75d0*Pi
c
c J_1(x) is odd in x, so we note the sign and use |x|
      sgn=1.d0
      if(x.lt.0.d0) sgn=-1.d0
      ax=dabs(x)
c
c use approximation appropriate to value of x:
c   range (0,4): rational approximation of J_1(x)
c   range [4,8): rational approximation of J_1(x)
c   range [8,12): rational approximation of
c                    J_1(x)-sqrt(2/(pi*x))*cos(x-3pi/4)
c   range [12,15): rational approximation of
c                    J_1(x)-sqrt(2/(pi*x))*cos(x-3pi/4)
c   range [15,20): rational approximation of
c                    J_1(x)-sqrt(2/(pi*x))*cos(x-3pi/4)
c   range [20,21): rational approximation of
c                    J_1(x)-sqrt(2/(pi*x))*cos(x-3pi/4)
c   range [21,infinity): Hankel's asymptotic expansion
c
      if (ax.eq.0.d0) then
        bj1=0.d0
      else if (ax.lt.4.d0) then
        call ratappr(ax,rabj1an,11,rabj1ad,11,bj1)
      else if (ax.lt.8.d0) then
        call ratappr(ax,rabj1bn,11,rabj1bd,11,bj1)
      else if (ax.lt.12.d0) then
        call ratappr(ax,rabj1cn,11,rabj1cd,11,bj1)
        bj1=bj1+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pi3d4)
      else if (ax.lt.15.d0) then
        call ratappr(ax,rabj1dn,11,rabj1dd,11,bj1)
        bj1=bj1+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pi3d4)
      else if (ax.lt.20.d0) then
        call ratappr(ax,rabj1en,11,rabj1ed,11,bj1)
        bj1=bj1+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pi3d4)
      else if (ax.lt.21.d0) then
        call ratappr(ax,rabj1fn,7,rabj1fd,7,bj1)
        bj1=bj1+dsqrt(2.d0/(Pi*ax))*dcos(ax-Pi3d4)
      else
        t=ax-Pi3d4
        call HankelPQ(1,ax,pp,qq)
        bj1=dsqrt(2.d0/(Pi*ax))*(pp*dcos(t)-qq*dsin(t))
      end if
      bj1=sgn*bj1
c
      return
      end
c
c***********************************************************************
      subroutine dbessjm(m,x,bjm)
c Bessel function of integer orders 0 through m-1: bjm(k) = J_{k-1}(x).
      implicit none
      integer m
      double precision x
      double precision bjm(*)
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      integer sgn,acc
      integer kk,ms,k,p
      double precision sd2,maxb,scl
      double precision bj0,bj1,bjp,bjc,bjn
      double precision ax,twovx,norm
c
c set sd2 to roughly sqr(# of desired significant digits)
c set maxb to renomalizing threshold
      sd2=169.d0
      maxb=1.d+16
      scl=1.d0/maxb
c
c sanity check
      if(m.le.0) then
        write(6,*) ' '
        write(6,*) ' <*** ERROR: dbessjm ***> requires argument m > 0.'
        write(6,*) '   Returning result array unchanged!'
        return
      end if
c
c note sgn(x) and use |x|; correct signs at end
      sgn=1
      if(x.lt.0.d0) sgn=-1
      ax=dabs(x)
      twovx=2.d0/ax
c
c clear result array
      do kk=1,m
        bjm(kk)=0.d0
      enddo
c
c treat x = 0 as a special case
      if(ax.eq.0.d0) then
        bjm(1)=1.d0
c use upward recursion if |x| exceeds maximum order
      else if(nint(ax).gt.m-1) then
        call dbessj0(x,bj0)
        bjm(1)=bj0
        if(m.eq.1) return
        call dbessj1(x,bj1)
        bjm(2)=bj1
        if(m.eq.2) return
        bjp=bj0
        bjc=bj1
        do k=1,m-2
          bjn=k*twovx*bjc-bjp
          bjp=bjc
          bjc=bjn
          bjm(k+2)=bjc
        end do
      else
c use downward recursion, starting well above desired maximum order
c   make ms even, set alternating flag acc to accumulate even J_k's,
c   and initialize variables
        ms=2*((m+nint(dsqrt(sd2*(m+1))))/2)
        acc=0
        norm=0.d0
        bjp=0.d0
        bjc=1.d0
        do k=ms,1,-1
          bjn=k*twovx*bjc-bjp
          bjp=bjc
          bjc=bjn
c   if results get out of hand, renormalize them
          if(dabs(bjc).gt.maxb) then
            norm=scl*norm
            bjp=scl*bjp
            bjc=scl*bjc
            do kk=1,m
              bjm(kk)=scl*bjm(kk)
            enddo
          end if
          if(k.le.m) bjm(k)=bjc
          if(acc.eq.1) then
            acc=0
            norm=norm+bjc
          else
            acc=1
          end if
        end do
c   recursion done; now normalize results using the identity
c     1 = J_0(x) + 2 J_2(x) + 2 J_4(x) + 2 J_6(x) + ...
        norm=1.d0/(2.d0*norm-bjc)
        do kk=1,m
          bjm(kk)=norm*bjm(kk)
        end do
      end if
c
c correct signs
      if(sgn.lt.0) then
        do kk=2,m,2
          bjm(kk)=-bjm(kk)
        end do
      end if
c
      return
      end
c
c***********************************************************************
      subroutine HankelPQ(m,x,P,Q)
c Compute the asymptotic functions P(m,x) and Q(m,x) that appear in
c Hankel's expansions of the Bessel functions for large x.
      implicit none
      integer m
      double precision x,P,Q
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      integer k,loopP,loopQ
      double precision mu,x8,e,sq,t,tp0,tp1,tq0,tq1
      double precision eps
c
c set tolerance:
c   this subroutine complains if the last term included in
c   the asymptotic series for P or Q exceeds this value
      eps=1.d-12
c
c initialize computation
      mu=4.d0*m*m
      x8=8.d0*x
      e=1.d0
      sq=e
      t=(mu-sq)/x8
      tp0=1.d0
      tq0=t
      P=tp0
      Q=tq0
      k=1
      loopP=1
      loopQ=1
c
c main loop
 1    continue
        tp1=-tp0*t
        k=k+1
        e=e+2
        sq=sq+e
        e=e+2
        sq=sq+e
        t=(mu-sq)/(k*x8)
        tp1=tp1*t
        tq1=tp1
        k=k+1
        e=e+2
        sq=sq+e
        e=e+2
        sq=sq+e
        t=(mu-sq)/(k*x8)
        tq1=tq1*t
        if (loopP.eq.1.and.dabs(tp1).lt.dabs(tp0)) then
          P=P+tp1
          tp0=tp1
        else
          if(loopP.eq.1.and.dabs(tp0).gt.eps) then
            write(6,10) eps
            write(6,11) tp0
          endif
          loopP=0
        end if
        if (loopQ.eq.1.and.dabs(tq1).lt.dabs(tq0)) then
          Q=Q+tq1
          tq0=tq1
        else
          if(loopQ.eq.1.and.dabs(tq0).gt.eps) then
            write(6,10) eps
            write(6,12) tq0
          endif
          loopQ=0
        end if
      if(loopP.eq.1.or.loopQ.eq.1) go to 1
c
c end main loop
c
 10   format('\n <*** WARNING: HankelPQ ***> tolerance ',1pe9.3)
 11   format(' not reached in HankelP series: last term = ',e10.3,'.')
 12   format(' not reached in HankelQ series: last term = ',e10.3,'.')
c
      return
      end
c
c***********************************************************************
      double precision function HankelP(m,x)
c Compute the asymptotic function P(m,x) that appears in Hankel's
c expansions of the Bessel functions for large x.
      implicit none
      integer m
      double precision x
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      integer k2,loop
      double precision mu,x82,t0,t1,e,sq,P
c
c initialize computation
      mu=4.d0*m*m
      x82=(8.d0*x)**2
      k2=2
      e=1
      sq=e
      t0=1.d0
      P=t0
      loop=1
c
c main loop
 1    continue
        t1=-t0/((k2-1)*k2*x82)
        t1=t1*(mu-sq)
        e=e+2
        sq=sq+e
        e=e+2
        sq=sq+e
        t1=t1*(mu-sq)
        if (dabs(t1).le.dabs(t0)) then
          P=P+t1
          t0=t1
          k2=k2+2
          e=e+2
          sq=sq+e
          e=e+2
          sq=sq+e
        else
          loop=0
        end if
      if(loop.eq.1) go to 1
c
c end main loop
c
      HankelP=P
      return
      end
c
c***********************************************************************
      double precision function HankelQ(m,x)
c Compute the asymptotic function Q(m,x) that appears in Hankel's
c expansions of the Bessel functions for large x.
      implicit none
      integer m
      double precision x
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      integer k2,loop
      double precision mu,x82,t0,t1,e,sq,Q
c
c initialize computation
      mu=4.d0*m*m
      x82=(8.d0*x)**2
      k2=2
      e=1
      sq=e
      t0=(mu-sq)/x82
      Q=t0
      loop=1
c
c main loop
 1    continue
        t1=-t0/(k2*(k2+1)*x82)
        e=e+2
        sq=sq+e
        e=e+2
        sq=sq+e
        t1=t1*(mu-sq)
        e=e+2
        sq=sq+e
        e=e+2
        sq=sq+e
        t1=t1*(mu-sq)
        if (dabs(t1).le.dabs(t0)) then
          Q=Q+t1
          t0=t1
          k2=k2+2
        else
          loop=0
        end if
      if(loop.eq.1) go to 1
c
c end main loop
c
      HankelQ=Q
      return
      end
c
c***********************************************************************
      subroutine ratappr(x,cfn,nn,cfd,nd,r)
c Compute the rational approximation for a function r(x)=rn(x)/rd(x),
c where the polynomials rn(x) and rd(x) are defined respectively by the
c nn coefficients in cfn and the nd coefficients in cfd.  The
c coefficients must be listed in order of increasing powers.  Thus, for
c the numerator, cfn(1) is the constant term, and cfn(nn) is the
c coefficient of the highest order term.  The polynomials are computed
c using Horner's method.  Return the result in r.
      implicit none
      integer nn,nd
      double precision x,r
      double precision cfn(*),cfd(*)
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      integer j
      double precision rn,rd
c
c compute numerator
      rn=cfn(nn)
      do j=nn-1,1,-1
        rn=cfn(j)+rn*x
      end do
c compute denominator
      rd=cfd(nd)
      do j=nn-1,1,-1
        rd=cfd(j)+rd*x
      end do
c compute ratio
      r=rn/rd
c
      return
      end
c
c***********************************************************************
