# GreenEOM

Ecuaciones de movimiento de funciones de Green a partir de un Hamiltoniano en
segunda cuantización, en Mathematica.

Le das H escrito con operadores de creación/aniquilación y una (o varias)
funciones de Green semilla. El paquete:

1. calcula `[A, H]` con las reglas fundamentales de (anti)conmutación,
2. escribe la ecuación de movimiento ya transformada a Fourier,
3. **cierra la jerarquía solo** — exacto si H es cuadrático,
4. resuelve el sistema algebraico lineal resultante,
5. construye `G^r`, `G^a` y `G^<` en equilibrio.

Es la generalización del esquema manual de la monografía de Judith (anillo
mesoscópico de Maiti): allí las ecuaciones se derivaban a mano y se tecleaban en
el notebook; aquí se derivan solas para cualquier H.

`GreenNEQ.wl` añade la capa de **no-equilibrio estacionario (Keldysh)**: contactos
continuos, `G^r`/`G^a`/`G^<` independientes, transmisión y corriente. Ver más
abajo.

## Por dónde empezar

Si llegas nuevo, en este orden:

1. **`examples.wl`** — nueve ejemplos trabajados, de un nivel resonante al
   límite atómico de Hubbard. Hace de manual: se lee mejor que esta sección.
2. **`tests.wl`** — ahí están las regresiones contra la monografía. El bloque
   `4a` reproduce la forma cerrada de tres sitios *carácter por carácter*,
   fases incluidas, comparándola contra la expresión transcrita en el propio
   script (`thesisG`, la ecuación numerada del artículo). Es lo que convence de
   que el generador no inventa, y no necesita ningún documento externo.
3. **`family.wl`** — los tres anillos (Judith, Emma, Dayanna) escritos como un
   solo Hamiltoniano con distintas filas encendidas.
4. **`open_family.wl`** — los mismos tres, pero con dos reservorios.

Todo corre en local:

```bash
cat GreenEOM.wl GreenNEQ.wl examples.wl > job.wl
wolframscript -file job.wl
```

Cada script imprime `PASS`/`FAIL` por línea y un recuento al final, así que si
algo se rompe se ve enseguida y se ve *dónde*.

## Estado

**Entorno de verificacion.** Mathematica 12.0 (kernel virgen por ssh a Astrum,
o wolframscript local); las comprobaciones de Python en `checks/` con Python 3
+ NumPy/SymPy (la simbolica exacta con SymPy 1.14, corrida tambien en Astrum).
Todo en doble precision; una cantidad se declara nula por debajo de 1e-10, y las
cifras reportadas como cero (1e-15 a 1e-11) quedan holgadamente bajo ese umbral.
Estados ocupados = los N autovalores mas bajos; los dos unicos puntos ambiguos
son los de nivel de Fermi casi degenerado, nombrados en el articulo.

**Verificar sin Wolfram.** El resultado cuantitativo central (la regla
N = 2 mod 4, la derivacion de la simetria, el barrido de 324 puntos) se
reproduce por completo en `checks/` con Python 3 + NumPy/SymPy/z3, sin
ninguna parte del paquete de Mathematica: `python checks/derivation.py`,
`python checks/exact_symmetry.py`, etc. El paquete `.wl` en si corre con
Wolfram Engine (gratuito para desarrolladores) via `wolframscript`, sin
necesidad de una licencia completa de Mathematica.

Verificado en kernel virgen (Mathematica 12), **122 pruebas**:

| Suite | Resultado |
|---|---|
| `tests.wl` (equilibrio) | **28/28** |
| `examples.wl` (ejemplos) | **9/9** |
| `tests_neq.wl` (no-equilibrio) | **24/24** |
| `transport.wl` (anillo abierto) | **16/16** |
| `impurities.wl` (robustez) | **13/13** |
| `family.wl` (grafeno + siliceno) | **13/13** |
| `open_family.wl` (los tres, abiertos) | **5/5** |
| `phsym.wl` (simetria particula-hueco) | **6/6** |
| `scans.wl` (barridos anchos) | **6/6** |
| `convention.wl` (convenios de angulo) | **2/2** |

Las pruebas de regresión más fuertes comparan la derivación simbólica contra
formas cerradas obtenidas a mano en la monografía y **transcritas literalmente en
`tests.wl`**: la expresión de referencia va impresa en el propio script, así que
la verificación es autocontenida y no requiere ningún documento externo (los
`.pdf`/`.nb` de la monografía son la procedencia de esas expresiones, no una
dependencia en tiempo de ejecución).

- **3 sitios**: reproduce carácter por carácter la forma cerrada de tres sitios
  (transcrita en `tests.wl` como `thesisG`, que es la ecuación numerada del
  artículo en §IV y procede de `tres sitios.pdf`), fases incluidas.
- **10 sitios**: reproduce `G_{n+1,n}(ω)` y `G_{n,n+1}(ω)` del capítulo de
  resultados (ecs. `green1` y siguiente).
- **2 sitios con Rashba + Dresselhaus**: la `eq1` generada coincide con la forma
  escrita a mano (transcrita en `tests.wl` como `want`, procedente de
  `TwoSites.nb`).
- **Límite atómico de Hubbard**: da el exacto `(1-n)/(ω-ε) + n/(ω-ε-U)`.
- **Nivel resonante**: produce la autoenergía de hibridación
  `Σ_k V_k²/(ω-ε_k)` sin que se le teclee.

## Uso

```mathematica
<< "GreenEOM.wl"

DeclareSpecies[c, "Fermion"];          (* por defecto ya es fermión *)
DeclareSpecies[b, "Boson"];

H = eps NCTimes[Cre[c, 1], Ann[c, 1]] + tt (NCTimes[Cre[c, 1], Ann[c, 2]] +
                                            NCTimes[Cre[c, 2], Ann[c, 1]]);

cl  = CloseEOM[{GF[Ann[c, 1], Cre[c, 1]]}, H];
sol = First[SolveGF[cl]];
G11 = GF[Mono[{Ann[c, 1]}], Mono[{Cre[c, 1]}]] /. sol
```

Los índices pueden ser cualquier expresión: `Ann[c, 3]`, `Ann[c, {n, s}]`,
`Ann[c, {sitio, espin, capa}]`. El delta de Kronecker sabe manejar índices
compuestos y desplazamientos simbólicos (`KroneckerDelta[n+1, n]` → 0).

### Piezas

| Función | Qué hace |
|---|---|
| `Ann[sp, i]`, `Cre[sp, i]` | operadores |
| `DeclareSpecies[sp, "Fermion"\|"Boson"]` | estadística (fermión por defecto) |
| `NCTimes[a, b, ...]` o `a ** b` | producto no conmutativo, lineal |
| `NOrder[e]` | orden normal, con sus contracciones |
| `Comm`, `AComm` | conmutador / anticonmutador |
| `GF[A, B]` | ⟨⟨A;B⟩⟩_ω, lineal en ambas ranuras |
| `EOMEquation[g, H]` | una sola ecuación de movimiento |
| `CloseEOM[semillas, H]` | cierra la jerarquía; devuelve ecuaciones + incógnitas |
| `SolveGF[cl]` | resuelve el sistema |
| `GFMatrix[cl]` | devuelve `{M, v}` con `M.G == v` — es `(ω − H)` |
| `SpectralDecomposition[g, ω]` | polos y residuos |
| `GLesser[g, ω, β, η]` | `G^<` en equilibrio |

### Opciones

| Opción | Por defecto | Para qué |
|---|---|---|
| `FrequencySymbol` | `ω` | símbolo de frecuencia |
| `SourceNormalization` | `1/Sqrt[2 Pi]` | convención de la transformada. Usa `1/(2 Pi)` para Zubarev, `1` para el resolvente desnudo |
| `SourceTime` | `0` | el `t'` de la fase `Exp[-I ω t']` |
| `MaxOperators` | `1` | orden máximo permitido. `1` = exacto para cuadráticos |
| `DecouplingRule` | `None` | función que trunca las GF que excedan `MaxOperators` |
| `MaxEquations` | `2000` | tope de seguridad contra índices infinitos |

## Lo que puede y lo que no

**Cierra exacto** cuando H es cuadrático: la EOM de una GF de una partícula solo
genera otras GF de una partícula, y un conjunto de índices finito las hace
finitas en número.

**No cierra** con interacción (Hubbard U, e-fonón, …): la jerarquía es infinita.
`CloseEOM` se **detiene con un error explícito** en vez de devolverte algo que
parece exacto y no lo es. Dos salidas:

- subir `MaxOperators` — a veces basta y sigue siendo exacto (el límite atómico
  de Hubbard cierra en 2 ecuaciones con `MaxOperators -> 3`);
- dar una `DecouplingRule` — tú eliges la aproximación y queda escrita en el
  código, no escondida.

**Requiere un conjunto de índices finito y explícito.** Construye H con
`Sum[..., {n, N}]` y `Mod` para las condiciones periódicas. Un H sobre índices
simbólicos no acotados dispara `CloseEOM::runaway`.

**Solo equilibrio / estacionario**, que es donde `w G = ... ` tiene sentido.
Fuera del equilibrio harían falta las dos ramas de Keldysh y las ecuaciones de
Dyson en el contorno; no está.

## No-equilibrio estacionario (`GreenNEQ.wl`)

### Qué cambia de verdad

La EOM **retardada** en estacionario tiene la misma estructura que en
equilibrio: $G^r=(\omega-H-\Sigma^r)^{-1}$. Lo que se rompe es otra cosa. En
equilibrio $G^<$ está atado a $G^r$ por un único factor de Fermi; **con dos
reservorios a potenciales distintos ese lazo se corta** y $G^<$ pasa a ser un
objeto independiente:

$$G^< = G^r\,\Sigma^<\,G^a,\qquad \Sigma^<=i\sum_\alpha \Gamma_\alpha f_\alpha$$

Y los contactos deben ser **continuos**: con niveles discretos solo hay polos,
nunca anchura finita ni corriente.

### Dos cosas se derivan, no se postulan

- **La autoenergía es el complemento de Schur.** `Downfold[h, keep, w]` elimina
  los grados de libertad fuera de `keep` y devuelve
  $\Sigma = H_{DL}(\omega-H_{LL})^{-1}H_{LD}$. Verificado contra la ruta EOM:
  da el mismo $\sum_k V_k^2/(\omega-\varepsilon_k)$.
- **Las reglas de Langreth son un motor de reescritura**, no una fórmula
  copiada. `ContourLesser[CProd[a,b,c]]` expande $(ABC)^<$ en sus tres términos,
  y aplicado a Dyson hace aparecer $G^r\Sigma^<G^a$.

### Uso

```mathematica
h    = HamiltonianMatrix[H, {Ann[c,1], Ann[c,2]}];   (* desde el conmutador *)
gamL = {{gL,0},{0,0}};  gamR = {{0,0},{0,gR}};
sig  = WideBandSelfEnergy[gamL + gamR];
gr   = RetardedG[h, sig, \[Omega]];
gl   = LesserG[gr, LesserSigma[{gamL,gamR}, {fL,fR}]];
T    = Transmission[gr, gamL, gamR];
```

### Validado contra

- **Nivel resonante**: $T(\omega)=\Gamma_L\Gamma_R/[(\omega-\varepsilon_0)^2+(\Gamma/2)^2]$
  exacto; $T=1$ en resonancia con acoplo simétrico; $T\le1$; regla de suma
  $\int A\,d\omega=1$.
- **El FDT se recupera, no se impone**: con $\mu_L=\mu_R$, $G^r\Sigma^<G^a$
  colapsa solo a $G^<=-(G^r-G^a)f$. Es la prueba de consistencia más fuerte de
  toda la capa.
- **Conservación de corriente** $I_L+I_R=0$ punto a punto en ω, y
  **Meir-Wingreen se reduce a Landauer** para H cuadrático.
- **Anillo AB con dos contactos**: reciprocidad de dos terminales
  $T(\phi)=T(-\phi)$, periodicidad en un cuanto de flujo, oscilación AB real.
- **Ocupación fuera de equilibrio**: simetría partícula-hueco
  $n(\mu)+n(-\mu)=1$, medio llenado con sesgo simétrico.

### Lo que NO cubre

Solo **H cuadrático, estacionario, sin correlaciones iniciales** — que es el
caso exactamente resoluble. No hay Keldysh interactuante (correcciones de
vértice), ni propagación temporal Kadanoff-Baym, ni régimen transitorio. Para
eso haría falta convolución en el contorno con dos tiempos, que es otro motor.

**Cuidado con la cola lorentziana** al integrar: $A(\omega)\sim\Gamma/2\pi\omega^2$
decae despacio, y una ventana de ±20 ya pierde ~3·10⁻³ del peso, suficiente para
romper las reglas de suma. Integra lejos y dale a `NIntegrate` los puntos de
resonancia.

## Resultado físico: el anillo de Judith, abierto

`transport.wl` conecta el anillo a dos contactos y calcula transporte
Landauer-Büttiker. `impurities.wl` intenta falsar el hallazgo.

Su monografía observa —sobre figuras de corrientes persistentes— que las
corrientes de espín Rashba y Dresselhaus son opuestas, y concluye que α=β debería
dar "una densidad muy próxima a cero". **En transporte eso es una identidad
exacta**, porque la permutación α↔β es una rotación de espín de π alrededor de
[110] que intercambia σx↔σy e invierte σz:

$$P_z(\alpha,\beta) = -P_z(\beta,\alpha)$$

de modo que α=β fuerza P_z=0 *exactamente*. Residuo peor 1.5·10⁻¹⁴ sobre 216
combinaciones y 3 tamaños de anillo.

**Y está protegida por simetría, no es ajuste fino.** Como U actúa solo en el
espacio de espín y es global, conmuta con cualquier término ∝1 en espín:

| Perturbación | Identidad |
|---|---|
| impureza de carga ε=0.9t | sobrevive (1e-15) |
| desorden aleatorio W=1t, W=10t | sobrevive (1e-15) |
| desorden + contactos fuera del eje | sobrevive (1e-15) |
| impureza Zeeman h=0.3t | **rompe** (0.77) |
| contactos polarizados Γ↑/Γ↓=2 | **rompe** (0.94) |
| Γ mayor pero igual en ambos espines | sobrevive (1e-15) |

Con σz no muere, **se deforma**: `P_z(α,β,h) = −P_z(β,α,−h)` se cumple a 2·10⁻¹⁵
donde la forma original falla por 0.77.

**Lectura experimental:** un P_z≠0 medido con acoplos balanceados es diagnóstico
directo de dispersión dependiente de espín — nada ciego al espín puede producirlo.

Ojo con los **dark states**: con contactos diametralmente opuestos los autoestados
impares respecto al eje tienen nodo en ambos contactos y `ω−H−Σ` es singular ahí.
La transmisión sí está bien definida; solo la inversa completa no. Un η pequeño lo
regulariza (hay test de independencia de η).

## Convenio del ángulo de enlace

Un salto de rango `m` que arranca en el sitio `n` lleva el punto medio del arco
que salva:

```
theta^(m)_n = phi_n + m Pi/N,      phi_n = 2 Pi (n-1)/N
```

Escrito así y no como el promedio `(phi_n + phi_{n+m})/2`, que parece lo mismo y
no lo es: en el enlace que cierra el anillo el promedio lee `phi_1` como `0` en
vez de `2 Pi` y deja ese enlace concreto a `Pi` de todos los demás. Con solo
espín-órbita de rango uno da igual — la diferencia es una constante, o sea una
rotación global de espín, invisible en cualquier observable polarizado en z —
pero deja de dar igual en cuanto hay términos de dos rangos distintos, que es
justo el anillo de siliceno. `convention.wl` mide las tres variantes lado a
lado: con el promedio ingenuo el resultado central del artículo desaparece.

**Excepción deliberada:** `ringH` en `tests.wl` mantiene el promedio de la
monografía, y lo dice en el propio código. Son regresiones contra funciones de
Green derivadas a mano en la tesis, y una regresión tiene que correr en el
convenio de aquello que reproduce; si no, prueba el modelo en vez del generador.

## Ejecutar las pruebas

Hay Mathematica en la máquina Windows, así que para cualquier cosa pequeña basta
concatenar y correr en local:

```bash
cat GreenEOM.wl GreenNEQ.wl tests.wl > job.wl
wolframscript -file job.wl
```

Para el cluster, `build_job.sh` concatena paquete + script en un solo fichero
(el bridge cuesta ~1 s por llamada, conviene agrupar) y redirige `Print` a un log
que se devuelve como valor:

```bash
bash build_job.sh tests.wl job.wl
```

Después, en un kernel virgen — que es la única verificación que cuenta, porque el
kernel del bridge es persistente y arrastra definiciones viejas:

```bash
scp job.wl astrum:~/greeneom_job.wl
ssh astrum '~/opt/Wolfram/Executables/wolframscript -file ~/greeneom_job.wl'
```

El fuente es **ASCII puro** (todo carácter especial va como `\[Nombre]`), para
que sobreviva el viaje PowerShell → ssh → shell remota sin mutilarse.

## Ficheros

- `GreenEOM.wl` — el paquete de equilibrio
- `GreenNEQ.wl` — capa de no-equilibrio estacionario (Keldysh)
- `tests.wl` — regresión de equilibrio (28)
- `examples.wl` — ejemplos trabajados, hace de manual (9)
- `tests_neq.wl` — regresión de no-equilibrio (24)
- `transport.wl` — anillo abierto con contactos, transporte Landauer (16)
- `impurities.wl` — falsación de la identidad P_z(α,β)=−P_z(β,α) (13)
- `family.wl` — anillos de grafeno (Emma Mora) y siliceno (Dayanna Pereira) (13)
- `open_family.wl` — los tres anillos abiertos con dos reservorios (5)
- `phsym.wl` — ruptura de la simetría partícula-hueco por los segundos vecinos
- `scans.wl` — los dos barridos anchos (324 contextos) de la Sec. III B
- `convention.wl` — compara los tres convenios de ángulo, lado a lado
- `figure4.wl`, `figures_family.wl` — generan los CSV de las figuras
- `checks/` — comprobaciones en Python, camino de codigo independiente del de
  Mathematica: `projective_algebra.py` (el lema {Pi, T^{N/2}} = 0 para
  N = 2 mod 4 en el anillo sin espin), `selection_rule.py` (si ese operador
  explica o no la regla en el anillo CON espin: no lo hace) y `sum_rule.py`
  (que es lo que de verdad se anula: la suma ocupada, no los niveles),
  `derivation.py` (la regla N = 2 mod 4 DERIVADA: traslacion generalizada
  T = T_1 x exp(-i sz pi/N), momentos semienteros, y las tres propiedades de
  los bloques; 7/7) y `exact_symmetry.py` (la conmutacion [T,H]=[T,V]=0 en
  SymPy EXACTO con acoplos simbolicos, no muestreo; corrida tambien en Astrum)
- `data/` — CSVs exportados y logs depositados
- `build_job.sh` — empaquetador para el bridge
