# Quimática - Orgánica: Documentación Completa

## Descripción General

Quimática - Orgánica es una aplicación web interactiva para la visualización y manipulación de estructuras moleculares. La aplicación permite a los usuarios visualizar moléculas en 2D y 3D, combinar estructuras moleculares usando IA, y exportar modelos 3D en varios formatos.

## Arquitectura del Sistema

### Backend (Flask)
- **Framework**: Flask 2.3.2
- **Procesamiento Molecular**: RDKit 2023.9.5
- **Visualización**: SVG para 2D, MOL blocks para 3D

### Frontend
- **Visualización 3D**: 3Dmol.js, Three.js
- **IA Generativa**: Google Generative AI (Gemini)
- **Interfaz**: HTML5, CSS3, JavaScript ES6 Modules

## Instalación y Configuración

### Requisitos del Sistema
- Python 3.8+
- Navegador web moderno con soporte para ES6 modules

### Instalación

```bash
# Clonar el repositorio
git clone <repository-url>
cd quimatica-organica

# Crear entorno virtual
python -m venv .venv
source .venv/bin/activate  # Linux/Mac
# o
.venv\Scripts\activate     # Windows

# Instalar dependencias
pip install -r requirements.txt

# Ejecutar la aplicación
python app.py
```

La aplicación estará disponible en `http://localhost:5000`

## API Backend

### Endpoints Principales

#### `GET /`
**Descripción**: Sirve la página principal de la aplicación.

**Respuesta**: Renderiza `templates/index.html`

**Ejemplo de uso**:
```bash
curl http://localhost:5000/
```

#### `GET /api/molecule/<molecule_name>`
**Descripción**: Obtiene la estructura de una molécula específica del banco de datos.

**Parámetros de consulta**:
- `show_atoms` (opcional): `true`/`false` - Mostrar índices de átomos en SVG (default: `true`)
- `show_bonds` (opcional): `true`/`false` - Mostrar índices de enlaces en SVG (default: `true`)

**Respuesta exitosa (200)**:
```json
{
  "mol": "string (MOL block format)",
  "smiles": "string (SMILES notation)",
  "svg": "string (SVG image)",
  "bonds": [
    {
      "bond_index": 0,
      "atom1_index": 1,
      "atom2_index": 2
    }
  ]
}
```

**Respuesta de error (404)**:
```json
{
  "error": "Molécula no encontrada"
}
```

**Ejemplo de uso**:
```bash
# Obtener estradiol con índices
curl "http://localhost:5000/api/molecule/estradiol?show_atoms=true&show_bonds=true"

# Obtener fulvestrant sin índices
curl "http://localhost:5000/api/molecule/fulvestrant?show_atoms=false&show_bonds=false"
```

#### `POST /api/render_smiles`
**Descripción**: Genera estructuras 2D y 3D a partir de una notación SMILES personalizada.

**Cuerpo de la petición**:
```json
{
  "smiles": "string (SMILES notation)",
  "show_atoms": boolean (opcional, default: true),
  "show_bonds": boolean (opcional, default: true)
}
```

**Respuesta exitosa (200)**:
```json
{
  "mol": "string (MOL block format)",
  "smiles": "string (SMILES notation)",
  "svg": "string (SVG image)",
  "bonds": [
    {
      "bond_index": 0,
      "atom1_index": 1,
      "atom2_index": 2
    }
  ]
}
```

**Respuesta de error (400)**:
```json
{
  "error": "No se proporcionó SMILES"
}
```

**Respuesta de error (422)**:
```json
{
  "error": "El SMILES proporcionado no es válido o no se pudo generar la estructura 3D."
}
```

**Ejemplo de uso**:
```bash
curl -X POST http://localhost:5000/api/render_smiles \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CCO",
    "show_atoms": true,
    "show_bonds": true
  }'
```

## Funciones Backend

### `smiles_to_mol_block(smiles)`
**Descripción**: Convierte una cadena SMILES a un bloque MOL 3D optimizado.

**Parámetros**:
- `smiles` (str): Notación SMILES de la molécula

**Retorna**:
- `str`: Bloque MOL 3D optimizado, o `None` si hay error

**Ejemplo**:
```python
mol_block = smiles_to_mol_block("CCO")  # Etanol
print(mol_block)  # Imprime el bloque MOL 3D
```

**Proceso interno**:
1. Convierte SMILES a objeto RDKit Mol
2. Añade hidrógenos explícitos
3. Genera conformación 3D usando ETKDG
4. Optimiza geometría con campo de fuerza MMFF94
5. Convierte a formato MOL block

### `smiles_to_svg(smiles, show_atom_indices=True, show_bond_indices=True)`
**Descripción**: Convierte SMILES a imagen SVG 2D y extrae información de enlaces.

**Parámetros**:
- `smiles` (str): Notación SMILES de la molécula
- `show_atom_indices` (bool): Mostrar índices de átomos
- `show_bond_indices` (bool): Mostrar índices de enlaces

**Retorna**:
- `tuple`: (svg_string, bonds_data) o (None, None) si hay error

**Ejemplo**:
```python
svg, bonds = smiles_to_svg("CCO", show_atom_indices=True, show_bond_indices=True)
print(f"SVG generado: {len(svg)} caracteres")
print(f"Enlaces encontrados: {len(bonds)}")
```

## Componentes Frontend

### Módulos JavaScript

#### `main.js` - Controlador Principal
**Descripción**: Punto de entrada de la aplicación, maneja la inicialización y eventos principales.

**Funciones principales**:

##### `init()`
**Descripción**: Inicializa todos los componentes de la aplicación.

**Proceso**:
1. Inicializa el visor 3D
2. Configura el menú contextual
3. Puebla las tarjetas de moléculas
4. Configura drag & drop
5. Registra event listeners

**Ejemplo de uso**:
```javascript
document.addEventListener('DOMContentLoaded', init);
```

##### `handleDownload()`
**Descripción**: Maneja la descarga de modelos 3D en diferentes formatos.

**Formatos soportados**:
- STL (Stereolithography)
- OBJ (Wavefront OBJ)
- FBX (Filmbox)
- 3MF (3D Manufacturing Format)
- AMF (Additive Manufacturing Format)

#### `api.js` - Cliente API
**Descripción**: Maneja la comunicación con servicios de IA generativa.

**Funciones exportadas**:

##### `getGenerativeModel(modelName)`
**Descripción**: Obtiene una instancia de modelo generativo de Google AI.

**Parámetros**:
- `modelName` (string): Nombre del modelo ('gemini-1.5-flash', 'gemini-1.5-pro', etc.)

**Retorna**:
- `GenerativeModel`: Instancia del modelo de IA

**Modelos disponibles**:
- `gemini-1.5-flash`: Rápido, para uso general
- `gemini-1.5-pro`: Más potente, para tareas complejas
- `gemini-2.5-flash-lite`: Ultrarrápido, para sugerencias
- `gemini-2.5-pro`: Última generación, máximo rendimiento

**Ejemplo**:
```javascript
import { getGenerativeModel } from './api.js';

const model = getGenerativeModel('gemini-1.5-flash');
const result = await model.generateContent('Describe esta molécula');
```

#### `viewer.js` - Visualización Molecular
**Descripción**: Maneja la visualización 2D y 3D de moléculas.

**Funciones principales**:

##### `initViewer()`
**Descripción**: Inicializa el visor 3D usando 3Dmol.js.

##### `loadMolecule(moleculeName)`
**Descripción**: Carga y visualiza una molécula por nombre.

**Parámetros**:
- `moleculeName` (string): Nombre de la molécula ('estradiol', 'fulvestrant')

##### `update3DStyle()`
**Descripción**: Actualiza el estilo de visualización 3D.

**Estilos disponibles**:
- `stick`: Representación de palos
- `line`: Estructura de alambre
- `sphere`: Esferas de van der Waals
- `ballandstick`: Combinación bolas y palos

##### `toggle3DLabels()`
**Descripción**: Alterna la visualización de etiquetas en el visor 3D.

#### `ui.js` - Interfaz de Usuario
**Descripción**: Maneja la interfaz de usuario y interacciones.

**Funciones principales**:

##### `populateMoleculeCards()`
**Descripción**: Puebla las tarjetas de moléculas en el banco molecular.

##### `updateCardSelection(moleculeName)`
**Descripción**: Actualiza la selección visual de tarjetas de moléculas.

##### `updateButtonState(text, loading, disabled)`
**Descripción**: Actualiza el estado de botones con indicadores de carga.

##### `initContextMenu()`
**Descripción**: Inicializa el menú contextual para sugerencias de IA.

##### `extractJson(rawText)`
**Descripción**: Extrae objetos JSON de respuestas de IA de forma robusta.

**Parámetros**:
- `rawText` (string): Respuesta cruda de la IA

**Retorna**:
- `object|null`: Objeto JSON parseado o null

#### `combination.js` - Combinación Molecular
**Descripción**: Maneja la lógica de combinación de moléculas usando IA.

**Funciones principales**:

##### `setupDragAndDrop()`
**Descripción**: Configura la funcionalidad de arrastrar y soltar para moléculas.

##### `combineMolecules()`
**Descripción**: Combina moléculas usando IA generativa.

##### `clearCombination(selectedMolecule)`
**Descripción**: Limpia la zona de combinación y restaura vista individual.

##### `handleSuggestionClick(suggestion)`
**Descripción**: Maneja clics en sugerencias de modificación molecular.

## Componentes de Interfaz

### Banco de Moléculas
**Descripción**: Sección que muestra las moléculas disponibles para visualización y combinación.

**Funcionalidades**:
- Visualización de tarjetas de moléculas
- Drag & drop para combinaciones
- Selección por clic para visualización individual

**Moléculas incluidas**:
- **Estradiol**: Hormona esteroide (SMILES: `C[C@]12CC[C@@H]3c4ccc(O)cc4CC[C@H]3[C@@H]1CC[C@@H]2O`)
- **Fulvestrant**: Antagonista del receptor de estrógeno

### Zona de Combinación
**Descripción**: Área interactiva para combinar moléculas usando IA.

**Componentes**:
- Dos slots para moléculas
- Selector de modelo de IA
- Botón de combinación/limpieza

**Modelos de IA disponibles**:
- Gemini 1.5 Flash (rápido)
- Gemini 1.5 Pro (avanzado)
- Gemini 2.5 Flash Lite (ultrarrápido)
- Gemini 2.5 Pro (última generación)

### Visor 2D
**Descripción**: Visualización bidimensional de estructuras moleculares.

**Controles**:
- Checkbox para índices de átomos
- Checkbox para índices de enlaces
- Pan y zoom interactivo

**Características**:
- Renderizado SVG de alta calidad
- Información de enlaces accesible
- Actualización en tiempo real

### Visor 3D
**Descripción**: Visualización tridimensional interactiva de moléculas.

**Controles disponibles**:
- **Estilos de representación**:
  - Palos (stick)
  - Alambre (line)
  - Esfera (sphere)
  - Bolas y palos (ballandstick)
- **Opciones de visualización**:
  - Mostrar/ocultar hidrógenos
  - Índices de átomos
  - Índices de enlaces
- **Controles de vista**:
  - Reiniciar vista
  - Rotación libre
  - Zoom

**Descarga de modelos 3D**:
- **Formatos soportados**: STL, OBJ, FBX, 3MF, AMF
- **Uso**: Ideal para impresión 3D y modelado

### Menú Contextual
**Descripción**: Menú emergente con sugerencias de IA para modificaciones moleculares.

**Funcionalidades**:
- Sugerencias automáticas basadas en contexto
- Input personalizado para modificaciones específicas
- Integración con modelos de IA
- Interfaz arrastrable

## Ejemplos de Uso

### Ejemplo 1: Visualización Básica de Molécula

```javascript
// Cargar y visualizar estradiol
import { loadMolecule } from './viewer.js';

// Cargar molécula
loadMolecule('estradiol');

// Cambiar estilo 3D
document.querySelector('input[value="sphere"]').checked = true;
update3DStyle();
```

### Ejemplo 2: Renderizado de SMILES Personalizado

```bash
# Renderizar cafeína
curl -X POST http://localhost:5000/api/render_smiles \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "show_atoms": true,
    "show_bonds": false
  }'
```

### Ejemplo 3: Combinación de Moléculas

```javascript
// Configurar combinación
const slot1 = document.getElementById('mol-slot-1');
const slot2 = document.getElementById('mol-slot-2');

// Simular drag & drop
slot1.dataset.molecule = 'estradiol';
slot2.dataset.molecule = 'fulvestrant';

// Ejecutar combinación
combineMolecules();
```

### Ejemplo 4: Descarga de Modelo 3D

```javascript
// Configurar formato de descarga
document.getElementById('download-format-select').value = 'stl';

// Descargar modelo
document.getElementById('download-btn').click();
```

## Personalización y Extensión

### Añadir Nuevas Moléculas

Para añadir moléculas al banco de datos:

```python
# En app.py
molecules_db = {
    "estradiol": "C[C@]12CC[C@@H]3c4ccc(O)cc4CC[C@H]3[C@@H]1CC[C@@H]2O",
    "fulvestrant": "C[C@]12CC[C@@H]3c4ccc(O)cc4CC[C@H]3[C@@H]1CC[C@H](C(F)(F)C(F)(F)S(=O)CCCCCCCCC)C2",
    "nueva_molecula": "SMILES_AQUI"  # Añadir nueva molécula
}
```

### Personalizar Estilos CSS

Los estilos se encuentran en `static/css/style.css`:

```css
/* Personalizar tarjetas de moléculas */
.molecule-card {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    border-radius: 15px;
    /* Más estilos... */
}
```

### Integrar Nuevos Modelos de IA

```javascript
// En api.js
export function getCustomModel(modelName) {
    return genAI.getGenerativeModel({ 
        model: modelName,
        generationConfig: {
            temperature: 0.7,
            topP: 0.8,
            maxOutputTokens: 1024,
        }
    });
}
```

## Solución de Problemas

### Problemas Comunes

#### Error: "Molécula no encontrada"
**Causa**: El nombre de molécula no existe en `molecules_db`
**Solución**: Verificar nombres disponibles o añadir nueva molécula

#### Error: "SMILES no válido"
**Causa**: Notación SMILES incorrecta o no reconocida por RDKit
**Solución**: Validar SMILES usando herramientas químicas

#### Visor 3D no carga
**Causa**: Problemas con 3Dmol.js o Three.js
**Solución**: 
1. Verificar conexión a internet
2. Comprobar consola del navegador
3. Actualizar navegador

#### IA no responde
**Causa**: Problemas con API key o límites de uso
**Solución**:
1. Verificar API key de Google AI
2. Comprobar cuotas de uso
3. Revisar conexión de red

### Logs y Debugging

Para habilitar logs detallados:

```python
# En app.py
import logging
logging.basicConfig(level=logging.DEBUG)

# Para debugging de RDKit
import rdkit.rdBase as rkrb
rkrb.LogToPythonLogger()
```

## Seguridad

### Consideraciones de Seguridad

1. **API Key**: La API key de Google AI está expuesta en el cliente. Para producción, implementar proxy backend.

2. **Validación de Input**: Validar todas las entradas SMILES para prevenir inyección de código.

3. **CORS**: Configurar CORS apropiadamente para producción.

4. **Rate Limiting**: Implementar límites de velocidad para las APIs.

### Recomendaciones para Producción

```python
# Configuración segura para producción
app.config['SECRET_KEY'] = 'your-secret-key-here'
app.config['DEBUG'] = False

# Implementar validación de SMILES
def validate_smiles(smiles):
    if len(smiles) > 1000:  # Límite de longitud
        return False
    # Más validaciones...
    return True
```

## Contribución

### Estructura de Archivos

```
quimatica-organica/
├── app.py                 # Aplicación Flask principal
├── requirements.txt       # Dependencias Python
├── templates/
│   └── index.html        # Template principal
├── static/
│   ├── css/
│   │   └── style.css     # Estilos principales
│   └── js/
│       ├── main.js       # Controlador principal
│       ├── api.js        # Cliente API
│       ├── viewer.js     # Visualización molecular
│       ├── ui.js         # Interfaz de usuario
│       └── combination.js # Combinación molecular
└── README.md             # Esta documentación
```

### Guías de Contribución

1. **Código Python**: Seguir PEP 8
2. **JavaScript**: Usar ES6+ y módulos
3. **CSS**: Usar metodología BEM
4. **Commits**: Mensajes descriptivos en español
5. **Testing**: Añadir tests para nuevas funcionalidades

### Roadmap

- [ ] Implementar más algoritmos de optimización molecular
- [ ] Añadir soporte para más formatos de archivo
- [ ] Integrar base de datos de moléculas más amplia
- [ ] Implementar colaboración en tiempo real
- [ ] Añadir análisis de propiedades moleculares
- [ ] Soporte para reacciones químicas

## Licencia

Este proyecto está bajo licencia MIT. Ver archivo LICENSE para más detalles.

## Contacto y Soporte

Para soporte técnico o preguntas sobre la aplicación:
- Crear issue en el repositorio
- Consultar documentación técnica
- Revisar ejemplos de código incluidos

---

*Documentación generada automáticamente para Quimática - Orgánica v1.0*