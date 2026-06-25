###############################################################
# TRABAJO DE FIN DE GRADO
# Relacion entre variables macroeconomicas espanolas (IPC y PIB)
# y el comportamiento bursatil de BBVA y Banco Santander (2005-2025)
#
# El script hace, en este orden:
#   PARTE 1 - Descarga y preparacion de los datos
#   PARTE 2 - Analisis exploratorio y diagnostico previo
#   PARTE 3 - Modelos: ARIMA/SARIMA, GARCH y cointegracion
#
# Las series financieras vienen de Yahoo Finance y las macro del INE.
###############################################################


###############################################################
# PARTE 1: DESCARGA Y PREPARACION DE DATOS
###############################################################

# Cargamos las librerias necesarias.
# Si alguna no esta instalada, instalarla una vez con:
#   install.packages(c("ineapir","quantmod","moments","FinTS","tseries",
#                      "forecast","rugarch","urca","vars","lmtest",
#                      "sandwich","xts","zoo"))
library(ineapir)    # descarga de datos del INE (IPC y PIB)
library(quantmod)   # descarga de precios de Yahoo Finance
library(moments)    # asimetria y curtosis
library(FinTS)      # test ARCH-LM
library(tseries)    # tests ADF, PP, KPSS, Jarque-Bera
library(forecast)   # auto.arima y previsiones
library(rugarch)    # modelos GARCH
library(urca)       # cointegracion (Engle-Granger, Phillips-Ouliaris, Johansen)
library(vars)       # seleccion de retardos (VARselect)
library(lmtest)     # test de Granger y contraste de coeficientes
library(sandwich)   # errores estandar Newey-West
library(xts)        # series temporales con fecha
library(zoo)        # apoyo a xts

# Guardamos todas las figuras en un unico PDF para la memoria.
pdf("figuras_TFG.pdf", width = 10, height = 7)
par(cex.main = 0.95)  # titulos un poco mas pequenos y consistentes en todo el PDF

# Las fechas del INE vienen en milisegundos; esta funcion las pasa a formato Date.
limpiar_fecha_ine <- function(x) {
  as.Date(as.POSIXct(as.numeric(x) / 1000, origin = "1970-01-01", tz = "UTC"))
}


##############################
# 1) IPC (INE) - DATOS MENSUALES
##############################

# La tabla 24077 del INE contiene el IPC general de Espana.
ipc_raw <- get_data_table(
  idTable   = 24077,
  dateStart = "2005/01/01",
  dateEnd   = "2025/12/31",
  unnest    = TRUE,
  validate  = FALSE
)

# La tabla puede traer varias series; nos quedamos con el IPC general nacional.
if ("Nombre" %in% names(ipc_raw)) {
  ipc_raw$Nombre <- trimws(as.character(ipc_raw$Nombre))
  filtro_general <- grepl("general", ipc_raw$Nombre, ignore.case = TRUE) &
    grepl("nacional|españa|espana", ipc_raw$Nombre, ignore.case = TRUE)
  if (sum(filtro_general, na.rm = TRUE) > 0) ipc_raw <- ipc_raw[filtro_general, ]
}

# Nos quedamos con fecha y valor, convertimos tipos y ordenamos.
ipc       <- ipc_raw[, c("Fecha", "Valor")]
ipc$Fecha <- limpiar_fecha_ine(ipc$Fecha)
ipc$Valor <- as.numeric(ipc$Valor)
ipc       <- ipc[order(ipc$Fecha), ]
ipc       <- ipc[ipc$Fecha >= as.Date("2005-01-01") & ipc$Fecha <= as.Date("2025-12-31"), ]
ipc       <- ipc[!duplicated(ipc$Fecha), ]   # por si quedara alguna fecha repetida

# Inflacion mensual = diferencia del logaritmo del IPC (tasa de variacion mensual).
ipc$infl_m <- c(NA, diff(log(ipc$Valor)))

# Inflacion interanual = IPC del mes respecto al mismo mes del ano anterior.
ipc$infl_yoy <- c(rep(NA, 12),
                  ipc$Valor[13:nrow(ipc)] / ipc$Valor[1:(nrow(ipc) - 12)] - 1)

cat("[IPC] Rango:", format(min(ipc$Fecha)), "a", format(max(ipc$Fecha)),
    "| Observaciones:", nrow(ipc), "\n")

# Pasamos las series a objeto ts (frecuencia mensual) para graficarlas.
ipc_ts       <- ts(ipc$Valor,    start = c(2005, 1), frequency = 12)
ipc_infl_m   <- ts(ipc$infl_m,   start = c(2005, 1), frequency = 12)
ipc_infl_yoy <- ts(ipc$infl_yoy, start = c(2005, 1), frequency = 12)

par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
plot(ipc_ts,       main = "IPC Espana (indice) 2005-2025",    ylab = "Indice",   xlab = "Ano")
plot(ipc_infl_m,   main = "Inflacion mensual (log-diff IPC)", ylab = "infl_m",   xlab = "Ano"); abline(h = 0)
plot(ipc_infl_yoy, main = "Inflacion interanual (YoY)",       ylab = "infl_yoy", xlab = "Ano"); abline(h = 0)
par(mfrow = c(1, 1))

# Base limpia del IPC que usaremos despues.
ipc_model <- ipc[, c("Fecha", "Valor", "infl_m", "infl_yoy")]


##############################
# 2) PIB (INE) - DATOS TRIMESTRALES
##############################

# La tabla 67822 contiene varias series de Contabilidad Nacional.
pib_raw <- get_data_table(
  idTable   = 67822,
  dateStart = "2005/01/01",
  dateEnd   = "2025/12/31",
  unnest    = TRUE,
  validate  = FALSE
)
pib_raw$Nombre <- trimws(as.character(pib_raw$Nombre))

# Elegimos el PIB a precios de mercado, en indices de volumen encadenados,
# dato base y ajustado de estacionalidad y calendario (asi medimos la economia real).
pib <- pib_raw[
  grepl("Datos ajustados de estacionalidad y calendario", pib_raw$Nombre, ignore.case = TRUE) &
  grepl("Producto interior bruto a precios de mercado",    pib_raw$Nombre, ignore.case = TRUE) &
  grepl("Dato base",                                       pib_raw$Nombre, ignore.case = TRUE) &
  grepl("Índices de volumen encadenados",                  pib_raw$Nombre, ignore.case = TRUE),
  c("Fecha", "Valor")
]

pib$Fecha <- limpiar_fecha_ine(pib$Fecha)
pib$Valor <- as.numeric(pib$Valor)
pib       <- pib[order(pib$Fecha), ]
pib       <- pib[pib$Fecha >= as.Date("2005-01-01") & pib$Fecha <= as.Date("2025-12-31"), ]
pib       <- pib[!duplicated(pib$Fecha), ]

# Crecimiento trimestral del PIB = diferencia del logaritmo (tasa de crecimiento).
pib$gdp_g_q <- c(NA, diff(log(pib$Valor)))

# Ano y trimestre de inicio para construir la serie ts trimestral.
start_year <- as.integer(format(min(pib$Fecha), "%Y"))
start_q    <- ((as.integer(format(min(pib$Fecha), "%m")) - 1) %/% 3) + 1

pib_ts   <- ts(pib$Valor,    start = c(start_year, start_q), frequency = 4)
pib_g_ts <- ts(pib$gdp_g_q,  start = c(start_year, start_q), frequency = 4)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(pib_ts,   main = "PIB Espana (trimestral) 2005-2025",         ylab = "Indice", xlab = "Ano")
plot(pib_g_ts, main = "Crecimiento trimestral del PIB (log-diff)", ylab = "crec_q", xlab = "Ano"); abline(h = 0)
par(mfrow = c(1, 1))


##############################
# 3) BANCOS: BBVA Y SANTANDER
##############################

# Esta funcion descarga un banco de Yahoo Finance y calcula precios,
# rendimientos y volatilidad. Se usa igual para BBVA y para Santander.
prep_bank <- function(ticker, nombre_banco) {

  # Descargamos los precios diarios. auto.assign = FALSE para que devuelva el objeto.
  x <- getSymbols(ticker, src = "yahoo",
                  from = "2005-01-01", to = "2026-01-01", auto.assign = FALSE)
  x <- na.omit(x)

  px_d     <- Ad(x)            # precio ajustado diario (incluye dividendos y splits)
  log_px_d <- log(px_d)        # log-precio diario
  px_m     <- to.monthly(px_d, indexAt = "lastof", OHLC = FALSE)  # precio mensual (ultimo del mes)
  log_px_m <- log(px_m)

  # Rendimientos logaritmicos: r_t = log(P_t) - log(P_{t-1}).
  r_d <- na.omit(dailyReturn(px_d, type = "log"))    # diarios
  r_m <- na.omit(monthlyReturn(px_m, type = "log"))  # mensuales

  # --- Graficos de niveles ---
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(index(px_d),     as.numeric(px_d),     type = "l", main = paste(nombre_banco, "- Precio diario"),     xlab = "Fecha", ylab = "Precio")
  plot(index(log_px_d), as.numeric(log_px_d), type = "l", main = paste(nombre_banco, "- Log(precio) diario"), xlab = "Fecha", ylab = "log(precio)")
  plot(index(px_m),     as.numeric(px_m),     type = "l", main = paste(nombre_banco, "- Precio mensual"),    xlab = "Fecha", ylab = "Precio")
  plot(index(log_px_m), as.numeric(log_px_m), type = "l", main = paste(nombre_banco, "- Log(precio) mensual"), xlab = "Fecha", ylab = "log(precio)")
  par(mfrow = c(1, 1))

  # --- Graficos de rendimientos ---
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(index(r_d), as.numeric(r_d), type = "l", main = paste(nombre_banco, "- Rendimiento diario"),  xlab = "Fecha", ylab = "r_d"); abline(h = 0)
  plot(index(r_m), as.numeric(r_m), type = "l", main = paste(nombre_banco, "- Rendimiento mensual"), xlab = "Fecha", ylab = "r_m"); abline(h = 0)
  hist(as.numeric(r_d), breaks = 60, main = paste(nombre_banco, "- Histograma rend. diarios"), xlab = "r_d")
  acf(as.numeric(r_d), main = paste(nombre_banco, "- ACF rendimientos diarios"))
  par(mfrow = c(1, 1))

  # --- Diagnostico de heterocedasticidad ---
  # Si la ACF de los rendimientos al cuadrado muestra autocorrelaciones, hay
  # efectos ARCH: la varianza no es constante. Eso justifica usar GARCH despues.
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf(as.numeric(r_d)^2, main = paste(nombre_banco, "- ACF de r_d^2"))
  plot(index(r_d), as.numeric(r_d), type = "l", main = paste(nombre_banco, "- Rendimientos diarios"), xlab = "Fecha", ylab = "r_d"); abline(h = 0)
  par(mfrow = c(1, 1))

  # Test ARCH-LM de Engle (1982). H0: no hay efectos ARCH.
  arch12 <- FinTS::ArchTest(as.numeric(r_d), lags = 12)
  arch24 <- FinTS::ArchTest(as.numeric(r_d), lags = 24)
  cat("\n[", nombre_banco, "] Test ARCH-LM (12 retardos):\n"); print(arch12)
  cat("\n[", nombre_banco, "] Test ARCH-LM (24 retardos):\n"); print(arch24)

  # --- Volatilidad realizada mensual ---
  # Desviacion tipica de los rendimientos diarios del mes, escalada por la raiz
  # del numero de sesiones del mes (asi pasa a volatilidad del periodo mensual).
  vol_m <- apply.monthly(r_d, function(z) sd(z, na.rm = TRUE) * sqrt(length(na.omit(z))))
  plot(index(vol_m), as.numeric(vol_m), type = "l",
       main = paste(nombre_banco, "- Volatilidad realizada mensual"), xlab = "Fecha", ylab = "vol_m")

  # Base mensual final del banco: rendimiento y volatilidad por mes.
  # Se unen por clave ano-mes porque el rendimiento se indexa al ultimo dia natural
  # y la volatilidad al ultimo dia con cotizacion.
  ret_m_df   <- data.frame(ym = format(as.Date(index(r_m)),   "%Y-%m"), ret_m = as.numeric(r_m))
  vol_m_df   <- data.frame(ym = format(as.Date(index(vol_m)), "%Y-%m"), vol_m = as.numeric(vol_m))
  monthly_df <- merge(ret_m_df, vol_m_df, by = "ym", all = FALSE)
  monthly_df$date <- as.Date(paste0(monthly_df$ym, "-01"))
  monthly_df <- monthly_df[order(monthly_df$date), c("date", "ret_m", "vol_m")]

  # Devolvemos los objetos para reutilizarlos en el resto del script.
  list(px_d = px_d, r_d = r_d, r_m = r_m, vol_m = vol_m, monthly_df = monthly_df)
}

# Ejecutamos la funcion para los dos bancos.
bbva <- prep_bank("BBVA.MC", "BBVA")
san  <- prep_bank("SAN.MC",  "SAN")


##############################
# 4) CONSTRUCCION DE LAS BASES FINALES
##############################

# Funciones auxiliares para crear claves de union temporal.
ym_from_date   <- function(d) format(as.Date(d), "%Y-%m")          # clave ano-mes "2020-03"
qkey_from_date <- function(d) {                                     # clave trimestre "2020Q1"
  y <- format(as.Date(d), "%Y")
  q <- ((as.integer(format(as.Date(d), "%m")) - 1) %/% 3) + 1
  paste0(y, "Q", q)
}

# --- 4.1 Bases mensuales: banco + IPC ---
ipc_model$ym <- ym_from_date(ipc_model$Fecha)

# BBVA + IPC
bbva_m       <- bbva$monthly_df
bbva_m$ym    <- ym_from_date(bbva_m$date)
data_m_bbva_ipc      <- merge(bbva_m, ipc_model[, c("ym", "infl_m", "infl_yoy")], by = "ym", all = FALSE)
data_m_bbva_ipc$date <- as.Date(paste0(data_m_bbva_ipc$ym, "-01"))
data_m_bbva_ipc      <- data_m_bbva_ipc[order(data_m_bbva_ipc$date), ]

# Santander + IPC
san_m        <- san$monthly_df
san_m$ym     <- ym_from_date(san_m$date)
data_m_san_ipc       <- merge(san_m, ipc_model[, c("ym", "infl_m", "infl_yoy")], by = "ym", all = FALSE)
data_m_san_ipc$date  <- as.Date(paste0(data_m_san_ipc$ym, "-01"))
data_m_san_ipc       <- data_m_san_ipc[order(data_m_san_ipc$date), ]

# --- 4.2 Bases trimestrales: banco + PIB ---
pib$qkey <- qkey_from_date(pib$Fecha)

# Rendimiento trimestral = suma de rendimientos diarios del trimestre (los log-rendimientos se suman).
# Volatilidad trimestral = desviacion tipica diaria escalada por la raiz de las sesiones del trimestre.
bbva_ret_q <- apply.quarterly(bbva$r_d, sum)
bbva_vol_q <- apply.quarterly(bbva$r_d, function(z) sd(z, na.rm = TRUE) * sqrt(length(na.omit(z))))
bbva_q     <- data.frame(date_q = as.Date(index(bbva_ret_q)),
                         ret_q = as.numeric(bbva_ret_q), vol_q = as.numeric(bbva_vol_q))
bbva_q$qkey     <- qkey_from_date(bbva_q$date_q)
data_q_bbva_pib <- merge(bbva_q, pib[, c("qkey", "gdp_g_q")], by = "qkey", all = FALSE)
data_q_bbva_pib <- data_q_bbva_pib[order(data_q_bbva_pib$date_q), ]

san_ret_q <- apply.quarterly(san$r_d, sum)
san_vol_q <- apply.quarterly(san$r_d, function(z) sd(z, na.rm = TRUE) * sqrt(length(na.omit(z))))
san_q     <- data.frame(date_q = as.Date(index(san_ret_q)),
                        ret_q = as.numeric(san_ret_q), vol_q = as.numeric(san_vol_q))
san_q$qkey     <- qkey_from_date(san_q$date_q)
data_q_san_pib <- merge(san_q, pib[, c("qkey", "gdp_g_q")], by = "qkey", all = FALSE)
data_q_san_pib <- data_q_san_pib[order(data_q_san_pib$date_q), ]

# --- 4.3 Bases mixtas: banco mensual + PIB trimestral ---
# Cada mes recibe el PIB del trimestre al que pertenece (el PIB se repite 3 veces).
# Solo se usan de forma descriptiva, no para modelizar.
bbva_m$qkey <- qkey_from_date(bbva_m$date)
data_mix_bbva_pib <- merge(bbva_m, pib[, c("qkey", "gdp_g_q")], by = "qkey", all = FALSE)
data_mix_bbva_pib <- data_mix_bbva_pib[order(data_mix_bbva_pib$date), ]

san_m$qkey <- qkey_from_date(san_m$date)
data_mix_san_pib  <- merge(san_m, pib[, c("qkey", "gdp_g_q")], by = "qkey", all = FALSE)
data_mix_san_pib  <- data_mix_san_pib[order(data_mix_san_pib$date), ]

cat("\nBases finales construidas:\n")
cat("  data_m_bbva_ipc   :", nrow(data_m_bbva_ipc),   "obs.\n")
cat("  data_m_san_ipc    :", nrow(data_m_san_ipc),    "obs.\n")
cat("  data_q_bbva_pib   :", nrow(data_q_bbva_pib),   "obs.\n")
cat("  data_q_san_pib    :", nrow(data_q_san_pib),    "obs.\n")
cat("  data_mix_bbva_pib :", nrow(data_mix_bbva_pib), "obs.\n")
cat("  data_mix_san_pib  :", nrow(data_mix_san_pib),  "obs.\n")

# --- 4.4 Variables dummy de crisis ---
# Valen 1 en los periodos de estres y 0 en el resto. Sirven para comprobar si
# la relacion entre macro y banca cambia en momentos excepcionales.
crear_dummies_crisis <- function(df, fecha_col) {
  f <- as.Date(df[[fecha_col]])
  df$crisis_2008         <- ifelse(f >= as.Date("2008-09-01") & f <= as.Date("2009-12-31"), 1, 0)
  df$covid_2020          <- ifelse(f >= as.Date("2020-03-01") & f <= as.Date("2020-12-31"), 1, 0)
  df$inflacion_2022_2023 <- ifelse(f >= as.Date("2022-01-01") & f <= as.Date("2023-12-31"), 1, 0)
  df
}
data_m_bbva_ipc <- crear_dummies_crisis(data_m_bbva_ipc, "date")
data_m_san_ipc  <- crear_dummies_crisis(data_m_san_ipc,  "date")
data_q_bbva_pib <- crear_dummies_crisis(data_q_bbva_pib, "date_q")
data_q_san_pib  <- crear_dummies_crisis(data_q_san_pib,  "date_q")


###############################################################
# PARTE 2: ANALISIS EXPLORATORIO Y DIAGNOSTICO PREVIO
###############################################################

##############################
# 5) REVISION GENERAL DE LAS BASES
##############################

# Muestra dimensiones, rango de fechas y valores perdidos de una base.
resumen_base <- function(df, nombre_base, fecha_var) {
  cat("\n==============================\n")
  cat("Base:", nombre_base, "\n")
  cat("Filas:", nrow(df), " | Columnas:", ncol(df), "\n")
  cat("Fecha minima:", as.character(min(df[[fecha_var]], na.rm = TRUE)), "\n")
  cat("Fecha maxima:", as.character(max(df[[fecha_var]], na.rm = TRUE)), "\n")
  cat("Valores perdidos por variable:\n"); print(colSums(is.na(df)))
}

resumen_base(data_m_bbva_ipc,   "MENSUAL_BBVA_IPC",     "date")
resumen_base(data_m_san_ipc,    "MENSUAL_SAN_IPC",      "date")
resumen_base(data_q_bbva_pib,   "TRIMESTRAL_BBVA_PIB",  "date_q")
resumen_base(data_q_san_pib,    "TRIMESTRAL_SAN_PIB",   "date_q")
resumen_base(data_mix_bbva_pib, "MIXTO_BBVA_PIB",       "date")
resumen_base(data_mix_san_pib,  "MIXTO_SAN_PIB",        "date")


##############################
# 6) GRAFICOS DE LAS SERIES FINALES
##############################

# BBVA + IPC
par(mfrow = c(4, 1), mar = c(4, 4, 3, 1))
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$ret_m,    type = "l", main = "BBVA - Rentabilidad mensual", xlab = "Fecha", ylab = "ret_m"); abline(h = 0)
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$vol_m,    type = "l", main = "BBVA - Volatilidad mensual",  xlab = "Fecha", ylab = "vol_m")
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$infl_m,   type = "l", main = "Inflacion mensual",           xlab = "Fecha", ylab = "infl_m"); abline(h = 0)
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$infl_yoy, type = "l", main = "Inflacion interanual",        xlab = "Fecha", ylab = "infl_yoy"); abline(h = 0)
par(mfrow = c(1, 1))

# SAN + IPC
par(mfrow = c(4, 1), mar = c(4, 4, 3, 1))
plot(data_m_san_ipc$date, data_m_san_ipc$ret_m,    type = "l", main = "SAN - Rentabilidad mensual", xlab = "Fecha", ylab = "ret_m"); abline(h = 0)
plot(data_m_san_ipc$date, data_m_san_ipc$vol_m,    type = "l", main = "SAN - Volatilidad mensual",  xlab = "Fecha", ylab = "vol_m")
plot(data_m_san_ipc$date, data_m_san_ipc$infl_m,   type = "l", main = "Inflacion mensual",          xlab = "Fecha", ylab = "infl_m"); abline(h = 0)
plot(data_m_san_ipc$date, data_m_san_ipc$infl_yoy, type = "l", main = "Inflacion interanual",       xlab = "Fecha", ylab = "infl_yoy"); abline(h = 0)
par(mfrow = c(1, 1))

# BBVA + PIB
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
plot(data_q_bbva_pib$date_q, data_q_bbva_pib$ret_q,   type = "l", main = "BBVA - Rentabilidad trimestral", xlab = "Fecha", ylab = "ret_q"); abline(h = 0)
plot(data_q_bbva_pib$date_q, data_q_bbva_pib$vol_q,   type = "l", main = "BBVA - Volatilidad trimestral",  xlab = "Fecha", ylab = "vol_q")
plot(data_q_bbva_pib$date_q, data_q_bbva_pib$gdp_g_q, type = "l", main = "Crecimiento trimestral del PIB", xlab = "Fecha", ylab = "gdp_g_q"); abline(h = 0)
par(mfrow = c(1, 1))

# SAN + PIB
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
plot(data_q_san_pib$date_q, data_q_san_pib$ret_q,   type = "l", main = "SAN - Rentabilidad trimestral", xlab = "Fecha", ylab = "ret_q"); abline(h = 0)
plot(data_q_san_pib$date_q, data_q_san_pib$vol_q,   type = "l", main = "SAN - Volatilidad trimestral",  xlab = "Fecha", ylab = "vol_q")
plot(data_q_san_pib$date_q, data_q_san_pib$gdp_g_q, type = "l", main = "Crecimiento trimestral del PIB", xlab = "Fecha", ylab = "gdp_g_q"); abline(h = 0)
par(mfrow = c(1, 1))

# Comparacion BBVA vs SAN (rentabilidad y volatilidad mensual)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$ret_m, type = "l", main = "Rentabilidad mensual: BBVA vs SAN", xlab = "Fecha", ylab = "ret_m")
lines(data_m_san_ipc$date, data_m_san_ipc$ret_m, lty = 2)
legend("topright", legend = c("BBVA", "SAN"), lty = c(1, 2), bty = "n"); abline(h = 0)
plot(data_m_bbva_ipc$date, data_m_bbva_ipc$vol_m, type = "l", main = "Volatilidad mensual: BBVA vs SAN", xlab = "Fecha", ylab = "vol_m")
lines(data_m_san_ipc$date, data_m_san_ipc$vol_m, lty = 2)
legend("topright", legend = c("BBVA", "SAN"), lty = c(1, 2), bty = "n")
par(mfrow = c(1, 1))


##############################
# 7) ESTADISTICOS DESCRIPTIVOS
##############################

# Estadisticos basicos de una variable (incluye asimetria y curtosis).
desc_stats <- function(x) {
  x <- x[!is.na(x)]
  c(n = length(x), media = mean(x), sd = sd(x), minimo = min(x), maximo = max(x),
    asimetria = moments::skewness(x), curtosis = moments::kurtosis(x))
}

desc_bbva_ipc <- rbind(ret_m = desc_stats(data_m_bbva_ipc$ret_m),  vol_m = desc_stats(data_m_bbva_ipc$vol_m),
                       infl_m = desc_stats(data_m_bbva_ipc$infl_m), infl_yoy = desc_stats(data_m_bbva_ipc$infl_yoy))
desc_san_ipc  <- rbind(ret_m = desc_stats(data_m_san_ipc$ret_m),   vol_m = desc_stats(data_m_san_ipc$vol_m),
                       infl_m = desc_stats(data_m_san_ipc$infl_m),  infl_yoy = desc_stats(data_m_san_ipc$infl_yoy))
desc_bbva_pib <- rbind(ret_q = desc_stats(data_q_bbva_pib$ret_q),  vol_q = desc_stats(data_q_bbva_pib$vol_q),
                       gdp_g_q = desc_stats(data_q_bbva_pib$gdp_g_q))
desc_san_pib  <- rbind(ret_q = desc_stats(data_q_san_pib$ret_q),   vol_q = desc_stats(data_q_san_pib$vol_q),
                       gdp_g_q = desc_stats(data_q_san_pib$gdp_g_q))

cat("\n--- Descriptivos BBVA + IPC ---\n"); print(round(desc_bbva_ipc, 6))
cat("\n--- Descriptivos SAN + IPC ---\n");  print(round(desc_san_ipc,  6))
cat("\n--- Descriptivos BBVA + PIB ---\n"); print(round(desc_bbva_pib, 6))
cat("\n--- Descriptivos SAN + PIB ---\n");  print(round(desc_san_pib,  6))


##############################
# 8) MATRICES DE CORRELACION
##############################

# Correlacion lineal contemporanea entre las variables de cada base.
cor_bbva_ipc <- cor(data_m_bbva_ipc[, c("ret_m", "vol_m", "infl_m", "infl_yoy")], use = "complete.obs")
cor_san_ipc  <- cor(data_m_san_ipc[,  c("ret_m", "vol_m", "infl_m", "infl_yoy")], use = "complete.obs")
cor_bbva_pib <- cor(data_q_bbva_pib[, c("ret_q", "vol_q", "gdp_g_q")], use = "complete.obs")
cor_san_pib  <- cor(data_q_san_pib[,  c("ret_q", "vol_q", "gdp_g_q")], use = "complete.obs")

cat("\n--- Correlaciones BBVA + IPC ---\n"); print(round(cor_bbva_ipc, 4))
cat("\n--- Correlaciones SAN + IPC ---\n");  print(round(cor_san_ipc,  4))
cat("\n--- Correlaciones BBVA + PIB ---\n"); print(round(cor_bbva_pib, 4))
cat("\n--- Correlaciones SAN + PIB ---\n");  print(round(cor_san_pib,  4))


##############################
# 9) DIAGNOSTICO PREVIO A LA MODELIZACION
##############################

# --- 9.1 ACF y PACF ---
# La ACF mide la correlacion con los retardos; la PACF descuenta los intermedios.
# Ayudan a identificar el orden de los modelos ARIMA.
graficos_acf_pacf <- function(x, nombre_variable) {
  x <- x[!is.na(x)]
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf(x,  main = paste("ACF -",  nombre_variable))
  pacf(x, main = paste("PACF -", nombre_variable))
  par(mfrow = c(1, 1))
}

graficos_acf_pacf(data_m_bbva_ipc$ret_m,    "BBVA ret_m")
graficos_acf_pacf(data_m_bbva_ipc$vol_m,    "BBVA vol_m")
graficos_acf_pacf(data_m_bbva_ipc$infl_m,   "IPC infl_m")
graficos_acf_pacf(data_m_bbva_ipc$infl_yoy, "IPC infl_yoy")
graficos_acf_pacf(data_q_bbva_pib$ret_q,    "BBVA ret_q")
graficos_acf_pacf(data_q_bbva_pib$vol_q,    "BBVA vol_q")
graficos_acf_pacf(data_q_bbva_pib$gdp_g_q,  "PIB gdp_g_q")

# --- 9.2 Tests de estacionariedad (ADF) y de normalidad (Jarque-Bera) ---
# ADF: H0 = raiz unitaria (serie no estacionaria).
# Jarque-Bera: H0 = los datos son normales.
tests_basicos <- function(x, nombre_variable) {
  x <- x[!is.na(x)]
  cat("\n==============================\n")
  cat("Variable:", nombre_variable, "\n")
  cat("Test ADF (H0: raiz unitaria):\n");        print(tseries::adf.test(x))
  cat("Test Jarque-Bera (H0: normalidad):\n");    print(tseries::jarque.bera.test(x))
}

tests_basicos(data_m_bbva_ipc$ret_m,    "BBVA ret_m")
tests_basicos(data_m_bbva_ipc$vol_m,    "BBVA vol_m")
tests_basicos(data_m_bbva_ipc$infl_m,   "IPC infl_m")
tests_basicos(data_m_bbva_ipc$infl_yoy, "IPC infl_yoy")
tests_basicos(data_q_bbva_pib$ret_q,    "BBVA ret_q")
tests_basicos(data_q_bbva_pib$vol_q,    "BBVA vol_q")
tests_basicos(data_q_bbva_pib$gdp_g_q,  "PIB gdp_g_q")
tests_basicos(data_m_san_ipc$ret_m,     "SAN ret_m")
tests_basicos(data_m_san_ipc$vol_m,     "SAN vol_m")
tests_basicos(data_q_san_pib$ret_q,     "SAN ret_q")
tests_basicos(data_q_san_pib$vol_q,     "SAN vol_q")


##############################
# 9.3 CAUSALIDAD DE GRANGER
##############################
# Comprueba si los valores pasados de una variable macro ayudan a predecir
# la rentabilidad o la volatilidad bancaria. OJO: es capacidad predictiva,
# no causalidad economica real.

# Test de Granger eligiendo el numero de retardos por AIC (con VARselect).
test_granger <- function(df, y_var, x_var, nombre, max_lag = 6) {
  aux <- na.omit(df[, c(y_var, x_var)]); names(aux) <- c("y", "x")
  if (nrow(aux) < (max_lag + 10)) {
    cat("\nGranger:", nombre, "- muestra demasiado pequena.\n")
    return(data.frame(contraste = nombre, y = y_var, x = x_var, lag = NA, F = NA, p_value = NA,
                      conclusion = "Muestra insuficiente", stringsAsFactors = FALSE))
  }
  lag_max_real <- max(1, min(max_lag, floor(nrow(aux) / 5)))
  sel <- vars::VARselect(aux, lag.max = lag_max_real, type = "const")
  p   <- as.integer(sel$selection["AIC(n)"]); if (is.na(p) || p < 1) p <- 1
  cat("\n==================================================\n")
  cat("Granger:", nombre, "| explica:", y_var, "| predictor:", x_var, "| retardos:", p, "\n")
  gt <- lmtest::grangertest(y ~ x, order = p, data = aux); print(gt)
  pval <- gt$`Pr(>F)`[2]; fval <- gt$F[2]
  conclusion <- ifelse(!is.na(pval) && pval < 0.05, "Se rechaza H0: hay evidencia predictiva al 5%", "No se rechaza H0 al 5%")
  data.frame(contraste = nombre, y = y_var, x = x_var, lag = p,
             F = round(fval, 5), p_value = round(pval, 5), conclusion = conclusion, stringsAsFactors = FALSE)
}

# Misma idea pero con un numero de retardos fijo (para ver si el resultado depende del lag).
test_granger_lag_fijo <- function(df, y_var, x_var, nombre, lag) {
  aux <- na.omit(df[, c(y_var, x_var)]); names(aux) <- c("y", "x")
  if (nrow(aux) < (lag + 10)) {
    return(data.frame(contraste = nombre, y = y_var, x = x_var, lag = lag, F = NA, p_value = NA,
                      conclusion = "Muestra insuficiente", stringsAsFactors = FALSE))
  }
  gt <- lmtest::grangertest(y ~ x, order = lag, data = aux)
  pval <- gt$`Pr(>F)`[2]; fval <- gt$F[2]
  conclusion <- ifelse(!is.na(pval) && pval < 0.05, "Se rechaza H0 al 5%", "No se rechaza H0 al 5%")
  data.frame(contraste = nombre, y = y_var, x = x_var, lag = lag,
             F = round(fval, 5), p_value = round(pval, 5), conclusion = conclusion, stringsAsFactors = FALSE)
}

granger_sensibilidad <- function(df, y_var, x_var, nombre, lags = 1:4) {
  do.call(rbind, lapply(lags, function(lag) test_granger_lag_fijo(df, y_var, x_var, nombre, lag)))
}

# Submuestra continua hasta un ano (para comprobar si el resultado depende del COVID).
submuestra_hasta_anio <- function(df, fecha_col, anio_final) {
  f <- as.Date(df[[fecha_col]])
  df[as.integer(format(f, "%Y")) <= anio_final, ]
}

cat("\n============================================================\n")
cat("RESULTADOS GRANGER\n")
cat("============================================================\n")

tabla_granger <- do.call(rbind, list(
  test_granger(data_m_bbva_ipc, "ret_m", "infl_m",  "IPC -> rentabilidad BBVA", max_lag = 6),
  test_granger(data_m_bbva_ipc, "vol_m", "infl_m",  "IPC -> volatilidad BBVA",  max_lag = 6),
  test_granger(data_m_san_ipc,  "ret_m", "infl_m",  "IPC -> rentabilidad SAN",  max_lag = 6),
  test_granger(data_m_san_ipc,  "vol_m", "infl_m",  "IPC -> volatilidad SAN",   max_lag = 6),
  test_granger(data_q_bbva_pib, "ret_q", "gdp_g_q", "PIB -> rentabilidad BBVA", max_lag = 4),
  test_granger(data_q_bbva_pib, "vol_q", "gdp_g_q", "PIB -> volatilidad BBVA",  max_lag = 4),
  test_granger(data_q_san_pib,  "ret_q", "gdp_g_q", "PIB -> rentabilidad SAN",  max_lag = 4),
  test_granger(data_q_san_pib,  "vol_q", "gdp_g_q", "PIB -> volatilidad SAN",   max_lag = 4)
))
cat("\nTABLA RESUMEN GRANGER:\n"); print(tabla_granger)

# Robustez 1: el PIB con varios retardos fijos.
cat("\n--- Sensibilidad Granger PIB con lags fijos ---\n")
tabla_granger_sensibilidad_pib <- do.call(rbind, list(
  granger_sensibilidad(data_q_bbva_pib, "ret_q", "gdp_g_q", "PIB -> rentabilidad BBVA", lags = 1:4),
  granger_sensibilidad(data_q_bbva_pib, "vol_q", "gdp_g_q", "PIB -> volatilidad BBVA",  lags = 1:4),
  granger_sensibilidad(data_q_san_pib,  "ret_q", "gdp_g_q", "PIB -> rentabilidad SAN",  lags = 1:4),
  granger_sensibilidad(data_q_san_pib,  "vol_q", "gdp_g_q", "PIB -> volatilidad SAN",   lags = 1:4)
))
print(tabla_granger_sensibilidad_pib)

# Robustez 1-bis: el IPC con varios retardos fijos (1 a 6).
# Comprueba que el resultado nulo del IPC no depende del numero de retardos.
# Genera la TABLA 5.3 (sensibilidad del IPC) de la memoria.
cat("\n--- Sensibilidad Granger IPC con lags fijos (1 a 6) ---\n")
tabla_granger_sensibilidad_ipc <- do.call(rbind, list(
  granger_sensibilidad(data_m_bbva_ipc, "ret_m", "infl_m", "IPC -> rentabilidad BBVA", lags = 1:6),
  granger_sensibilidad(data_m_bbva_ipc, "vol_m", "infl_m", "IPC -> volatilidad BBVA",  lags = 1:6),
  granger_sensibilidad(data_m_san_ipc,  "ret_m", "infl_m", "IPC -> rentabilidad SAN",  lags = 1:6),
  granger_sensibilidad(data_m_san_ipc,  "vol_m", "infl_m", "IPC -> volatilidad SAN",   lags = 1:6)
))
print(tabla_granger_sensibilidad_ipc)

# Robustez 2: el PIB usando solo datos hasta 2019 (sin el shock COVID).
cat("\n--- Granger PIB en submuestra hasta 2019 ---\n")
data_q_bbva_pib_pre2020 <- submuestra_hasta_anio(data_q_bbva_pib, "date_q", 2019)
data_q_san_pib_pre2020  <- submuestra_hasta_anio(data_q_san_pib,  "date_q", 2019)
tabla_granger_pib_pre2020 <- do.call(rbind, list(
  test_granger(data_q_bbva_pib_pre2020, "ret_q", "gdp_g_q", "PIB -> rentabilidad BBVA hasta 2019", max_lag = 4),
  test_granger(data_q_bbva_pib_pre2020, "vol_q", "gdp_g_q", "PIB -> volatilidad BBVA hasta 2019",  max_lag = 4),
  test_granger(data_q_san_pib_pre2020,  "ret_q", "gdp_g_q", "PIB -> rentabilidad SAN hasta 2019",  max_lag = 4),
  test_granger(data_q_san_pib_pre2020,  "vol_q", "gdp_g_q", "PIB -> volatilidad SAN hasta 2019",   max_lag = 4)
))
print(tabla_granger_pib_pre2020)


##############################
# 9.4 REGRESIONES CON DUMMIES DE CRISIS
##############################
# Regresamos la variable bancaria sobre la macro, las dummies de crisis y sus
# interacciones. Las interacciones dicen si la pendiente cambia en crisis.
# Usamos errores Newey-West, robustos a autocorrelacion y heterocedasticidad.

ajustar_modelo_crisis <- function(df, y_var, x_var, nombre) {
  aux <- na.omit(df[, c(y_var, x_var, "crisis_2008", "covid_2020", "inflacion_2022_2023")])
  formula_txt <- paste0(y_var, " ~ ", x_var,
                        " + crisis_2008 + covid_2020 + inflacion_2022_2023",
                        " + ", x_var, ":crisis_2008 + ", x_var, ":covid_2020 + ", x_var, ":inflacion_2022_2023")
  cat("\n==================================================\n")
  cat("Modelo con dummies de crisis:", nombre, "\n")
  mod <- lm(as.formula(formula_txt), data = aux)
  print(summary(mod))
  # Numero de retardos de Newey-West segun una regla habitual que depende de n.
  lag_nw <- max(1, floor(4 * (nobs(mod) / 100)^(2 / 9)))
  rob <- lmtest::coeftest(mod, vcov. = sandwich::NeweyWest(mod, lag = lag_nw, prewhite = FALSE))
  cat("\nCoeficientes con errores Newey-West (lag =", lag_nw, "):\n"); print(rob)
  if (nobs(mod) < 100) cat("Nota: muestra corta para tantas dummies e interacciones; interpretar con cautela.\n")
  data.frame(modelo = nombre, n = nobs(mod), lag_Newey_West = lag_nw,
             R2_ajustado = round(summary(mod)$adj.r.squared, 5), stringsAsFactors = FALSE)
}

tabla_crisis <- do.call(rbind, list(
  ajustar_modelo_crisis(data_m_bbva_ipc, "ret_m", "infl_m",  "BBVA rentabilidad mensual + IPC"),
  ajustar_modelo_crisis(data_m_bbva_ipc, "vol_m", "infl_m",  "BBVA volatilidad mensual + IPC"),
  ajustar_modelo_crisis(data_m_san_ipc,  "ret_m", "infl_m",  "SAN rentabilidad mensual + IPC"),
  ajustar_modelo_crisis(data_m_san_ipc,  "vol_m", "infl_m",  "SAN volatilidad mensual + IPC"),
  ajustar_modelo_crisis(data_q_bbva_pib, "ret_q", "gdp_g_q", "BBVA rentabilidad trimestral + PIB"),
  ajustar_modelo_crisis(data_q_bbva_pib, "vol_q", "gdp_g_q", "BBVA volatilidad trimestral + PIB"),
  ajustar_modelo_crisis(data_q_san_pib,  "ret_q", "gdp_g_q", "SAN rentabilidad trimestral + PIB"),
  ajustar_modelo_crisis(data_q_san_pib,  "vol_q", "gdp_g_q", "SAN volatilidad trimestral + PIB")
))
cat("\nTABLA RESUMEN REGRESIONES CON CRISIS:\n"); print(tabla_crisis)


###############################################################
# PARTE 3: MODELOS (ARIMA/SARIMA, GARCH Y COINTEGRACION)
###############################################################

##############################
# FUNCIONES AUXILIARES
##############################

# Pasa una clave trimestral "2020Q1" a la fecha del primer dia del trimestre.
qkey_to_date <- function(qkey) {
  y <- as.integer(substr(qkey, 1, 4)); q <- as.integer(substr(qkey, 6, 6))
  as.Date(paste0(y, "-", sprintf("%02d", (q - 1) * 3 + 1), "-01"))
}

# Saca el p-valor de un test (devuelve NA si no hay).
extraer_p <- function(test_obj) if (inherits(test_obj, "htest")) test_obj$p.value else NA_real_

# Ejecuta algo y, si da error, devuelve NULL en vez de parar el script.
safe_test <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)

# Test de Ljung-Box (H0: no hay autocorrelacion en los residuos) con proteccion
# para muestras cortas. fitdf descuenta los parametros AR y MA del modelo.
box_ljung_seguro <- function(res, lag_max = 12, fitdf = 0) {
  res <- na.omit(as.numeric(res)); n <- length(res)
  if (n <= 10) return(NULL)
  lag_usado <- min(lag_max, n - 1)
  if (lag_usado <= fitdf) lag_usado <- min(n - 1, fitdf + 1)
  if (lag_usado <= 0 || lag_usado >= n) return(NULL)
  Box.test(res, lag = lag_usado, type = "Ljung-Box", fitdf = fitdf)
}

# Crea un objeto ts a partir de un data.frame y una columna de fecha.
crear_ts <- function(df, fecha_col, var_col, frecuencia) {
  aux <- df[, c(fecha_col, var_col)]; names(aux) <- c("fecha", "valor")
  aux$fecha <- as.Date(aux$fecha); aux$valor <- as.numeric(aux$valor)
  aux <- aux[order(aux$fecha), ]; aux <- aux[!is.na(aux$valor), ]
  start_year <- as.integer(format(min(aux$fecha), "%Y"))
  start_per  <- if (frecuencia == 12) as.integer(format(min(aux$fecha), "%m"))
                else ((as.integer(format(min(aux$fecha), "%m")) - 1) %/% 3) + 1
  ts(aux$valor, start = c(start_year, start_per), frequency = frecuencia)
}


###############################################################
# 10) MODELOS ARIMA / SARIMA
###############################################################
# Modelizamos la dinamica EN MEDIA de las series ya transformadas
# (rendimientos, volatilidades, inflacion, crecimiento del PIB).
# A los precios en nivel NO se les aplica ARIMA: eso se deja para cointegracion.

# --- 10.1 Tests de raiz unitaria previos (ADF, PP y KPSS) ---
# Los tres juntos dan una clasificacion mas prudente.
# ADF y PP: H0 = raiz unitaria. KPSS: H0 = estacionariedad (justo al reves).
test_raiz_unitaria <- function(x, nombre) {
  x <- na.omit(as.numeric(x))
  cat("\n--------------------------------------------------\n")
  cat("Tests de raiz unitaria:", nombre, "\n")
  adf  <- safe_test(tseries::adf.test(x))
  pp   <- safe_test(tseries::pp.test(x))
  kpss <- safe_test(tseries::kpss.test(x, null = "Level"))
  adf_p <- extraer_p(adf); pp_p <- extraer_p(pp); kpss_p <- extraer_p(kpss)
  cat("ADF  p-value:", round(adf_p,  4), " | H0: raiz unitaria\n")
  cat("PP   p-value:", round(pp_p,   4), " | H0: raiz unitaria\n")
  cat("KPSS p-value:", round(kpss_p, 4), " | H0: estacionariedad\n")
  conclusion <- ifelse(!is.na(adf_p) && adf_p < 0.05,
                       "ADF rechaza raiz unitaria: serie estacionaria",
                       "ADF no rechaza raiz unitaria: revisar diferenciacion")
  cat("Conclusion:", conclusion, "\n")
  data.frame(serie = nombre, ADF_p = round(adf_p, 5), PP_p = round(pp_p, 5),
             KPSS_p = round(kpss_p, 5), conclusion = conclusion, stringsAsFactors = FALSE)
}

# --- 10.2 Ajuste ARIMA con diagnostico de residuos y prevision ---
# auto.arima elige el orden (p,d,q) minimizando el AIC.
# Para rendimientos, volatilidades y PIB usamos ARIMA NO estacional.
ajustar_arima <- function(x_ts, nombre, h) {
  x_ts <- na.omit(x_ts)
  cat("\n==================================================\n")
  cat("ARIMA para:", nombre, "\n")
  modelo <- forecast::auto.arima(x_ts, seasonal = FALSE, stepwise = FALSE,
                                 approximation = FALSE, ic = "aic", trace = FALSE)
  print(summary(modelo))
  res   <- na.omit(as.numeric(residuals(modelo)))
  orden <- forecast::arimaorder(modelo)
  fitdf <- orden[1] + orden[3]               # p + q (parametros AR y MA)
  lb12  <- box_ljung_seguro(res, 12, fitdf)
  lb24  <- box_ljung_seguro(res, 24, fitdf)
  jb    <- safe_test(tseries::jarque.bera.test(res))
  lb12_p <- extraer_p(lb12); lb24_p <- extraer_p(lb24); jb_p <- extraer_p(jb)
  cat("\n--- Diagnostico de residuos ---\n")
  cat("Ljung-Box lag 12 p-value:", round(lb12_p, 4), " (H0: sin autocorrelacion)\n")
  cat("Ljung-Box lag 24 p-value:", round(lb24_p, 4), "\n")
  cat("Jarque-Bera   p-value   :", round(jb_p,   4), " (H0: normalidad)\n")
  cat("AIC:", round(AIC(modelo), 4), " | BIC:", round(BIC(modelo), 4), "\n")

  # Graficos de diagnostico de los residuos.
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(res, type = "l", main = paste("Residuos ARIMA -", nombre), xlab = "Observacion", ylab = "residuo"); abline(h = 0)
  acf(res, main = paste("ACF residuos -", nombre))
  hist(res, breaks = 30, probability = TRUE, main = paste("Histograma residuos -", nombre), xlab = "residuo")
  curve(dnorm(x, mean = mean(res), sd = sd(res)), add = TRUE, lwd = 2)
  qqnorm(res, main = paste("QQ-plot residuos -", nombre)); qqline(res, lwd = 2)
  par(mfrow = c(1, 1))

  # Prevision ilustrativa (no es validacion fuera de muestra).
  prev <- forecast::forecast(modelo, h = h)
  plot(prev, main = paste("Prevision ARIMA -", nombre), xlab = "Periodo", ylab = "valor")

  if ((!is.na(lb12_p) && lb12_p < 0.05) || (!is.na(lb24_p) && lb24_p < 0.05))
    cat("AVISO: el ARIMA deja autocorrelacion residual significativa.\n")

  resumen <- data.frame(
    serie = nombre, modelo = paste0("ARIMA(", orden[1], ",", orden[2], ",", orden[3], ")"),
    AIC = round(AIC(modelo), 4), BIC = round(BIC(modelo), 4),
    LjungBox_p12 = round(lb12_p, 5), LjungBox_p24 = round(lb24_p, 5),
    JarqueBera_p = round(jb_p, 5), horizonte_prevision = h, stringsAsFactors = FALSE)
  list(modelo = modelo, residuos = res, forecast = prev, resumen = resumen)
}

# --- 10.3 Series que vamos a modelizar ---
series_tests_raiz <- list(
  list(nombre = "BBVA_ret_m",   x = crear_ts(data_m_bbva_ipc, "date",   "ret_m",    12), h = 12),
  list(nombre = "SAN_ret_m",    x = crear_ts(data_m_san_ipc,  "date",   "ret_m",    12), h = 12),
  list(nombre = "BBVA_vol_m",   x = crear_ts(data_m_bbva_ipc, "date",   "vol_m",    12), h = 12),
  list(nombre = "SAN_vol_m",    x = crear_ts(data_m_san_ipc,  "date",   "vol_m",    12), h = 12),
  list(nombre = "IPC_infl_m",   x = crear_ts(data_m_bbva_ipc, "date",   "infl_m",   12), h = 12),
  list(nombre = "IPC_infl_yoy", x = crear_ts(data_m_bbva_ipc, "date",   "infl_yoy", 12), h = 12),
  list(nombre = "BBVA_ret_q",   x = crear_ts(data_q_bbva_pib, "date_q", "ret_q",     4), h =  4),
  list(nombre = "SAN_ret_q",    x = crear_ts(data_q_san_pib,  "date_q", "ret_q",     4), h =  4),
  list(nombre = "PIB_gdp_g_q",  x = crear_ts(data_q_bbva_pib, "date_q", "gdp_g_q",   4), h =  4)
)

# Para la inflacion NO ponemos ARIMA no estacional en la tabla final: como el IPC
# mensual puede tener estacionalidad, ahi usamos SARIMA (mas abajo).
series_arima <- list(
  list(nombre = "BBVA_ret_m",  x = crear_ts(data_m_bbva_ipc, "date",   "ret_m",  12), h = 12),
  list(nombre = "SAN_ret_m",   x = crear_ts(data_m_san_ipc,  "date",   "ret_m",  12), h = 12),
  list(nombre = "BBVA_vol_m",  x = crear_ts(data_m_bbva_ipc, "date",   "vol_m",  12), h = 12),
  list(nombre = "SAN_vol_m",   x = crear_ts(data_m_san_ipc,  "date",   "vol_m",  12), h = 12),
  list(nombre = "BBVA_ret_q",  x = crear_ts(data_q_bbva_pib, "date_q", "ret_q",   4), h =  4),
  list(nombre = "SAN_ret_q",   x = crear_ts(data_q_san_pib,  "date_q", "ret_q",   4), h =  4),
  list(nombre = "PIB_gdp_g_q", x = crear_ts(data_q_bbva_pib, "date_q", "gdp_g_q", 4), h =  4)
)

# --- 10.4 Ejecucion ---
cat("\n============================================================\n")
cat("RESULTADOS ARIMA\n")
cat("============================================================\n")

cat("\n=== TESTS DE RAIZ UNITARIA ===\n")
tabla_raiz_unitaria <- do.call(rbind, lapply(series_tests_raiz, function(s) test_raiz_unitaria(s$x, s$nombre)))

cat("\n=== MODELOS ARIMA ===\n")
resultados_arima <- lapply(series_arima, function(s) ajustar_arima(s$x, s$nombre, s$h))
tabla_arima_no_inflacion <- do.call(rbind, lapply(resultados_arima, function(r) r$resumen))
cat("\nTABLA ARIMA NO ESTACIONAL (sin inflacion):\n"); print(tabla_arima_no_inflacion)

# --- 10.5 SARIMA para la inflacion ---
# auto.arima con seasonal = TRUE permite parte estacional. La inflacion mensual e
# interanual encajan mejor aqui: bajan el AIC y dejan residuos sin autocorrelacion.
ajustar_sarima_inflacion <- function(x_ts, nombre, h = 12) {
  x_ts <- na.omit(x_ts)
  cat("\n==================================================\n")
  cat("SARIMA para inflacion:", nombre, "\n")
  modelo <- forecast::auto.arima(x_ts, seasonal = TRUE, stepwise = FALSE,
                                 approximation = FALSE, ic = "aic", trace = FALSE)
  print(summary(modelo))
  res   <- na.omit(as.numeric(residuals(modelo)))
  orden <- forecast::arimaorder(modelo)
  fitdf <- sum(orden[c("p", "q", "P", "Q")], na.rm = TRUE)
  lb12  <- box_ljung_seguro(res, 12, fitdf); lb24 <- box_ljung_seguro(res, 24, fitdf)
  jb    <- safe_test(tseries::jarque.bera.test(res))
  lb12_p <- extraer_p(lb12); lb24_p <- extraer_p(lb24); jb_p <- extraer_p(jb)
  cat("Ljung-Box lag 12 p-value:", round(lb12_p, 4), "\n")
  cat("Ljung-Box lag 24 p-value:", round(lb24_p, 4), "\n")
  cat("Jarque-Bera   p-value   :", round(jb_p,   4), "\n")
  frec <- if ("Frequency" %in% names(orden)) orden["Frequency"] else frequency(x_ts)
  data.frame(
    serie = sub("_SARIMA$", "", nombre),
    modelo = paste0("SARIMA ARIMA(", orden["p"], ",", orden["d"], ",", orden["q"], ")(",
                    orden["P"], ",", orden["D"], ",", orden["Q"], ")[", frec, "]"),
    AIC = round(AIC(modelo), 4), BIC = round(BIC(modelo), 4),
    LjungBox_p12 = round(lb12_p, 5), LjungBox_p24 = round(lb24_p, 5),
    JarqueBera_p = round(jb_p, 5), horizonte_prevision = h, stringsAsFactors = FALSE)
}

cat("\n=== SARIMA PARA INFLACION ===\n")
tabla_sarima_inflacion <- rbind(
  ajustar_sarima_inflacion(crear_ts(data_m_bbva_ipc, "date", "infl_m",   12), "IPC_infl_m_SARIMA",   12),
  ajustar_sarima_inflacion(crear_ts(data_m_bbva_ipc, "date", "infl_yoy", 12), "IPC_infl_yoy_SARIMA", 12)
)
print(tabla_sarima_inflacion)

# Tabla ARIMA/SARIMA final: ARIMA no estacional para todo menos inflacion, y SARIMA para inflacion.
tabla_arima <- rbind(tabla_arima_no_inflacion, tabla_sarima_inflacion)
rownames(tabla_arima) <- NULL
cat("\nTABLA PRINCIPAL ARIMA/SARIMA:\n"); print(tabla_arima)


###############################################################
# 11) MODELOS GARCH
###############################################################
# Modelizamos la volatilidad condicional DIARIA de cada banco.
# Comparamos cuatro modelos y nos quedamos con el mejor segun BIC:
#   1) GARCH(1,1) normal
#   2) GARCH(1,1) t de Student            -> colas pesadas
#   3) GJR-GARCH(1,1) t                   -> asimetria (efecto apalancamiento)
#   4) EGARCH(1,1) t                      -> asimetria sobre el log de la varianza

# --- 11.1 Preparacion de los rendimientos diarios ---
# Reutilizamos los rendimientos diarios de la Parte 1 y los pasamos a porcentaje
# (multiplicar por 100 mejora el comportamiento numerico del optimizador).
preparar_rend_diario <- function(r_d_xts, nombre_banco) {
  rd_xts <- na.omit(r_d_xts) * 100
  colnames(rd_xts) <- "ret_d"
  cat("[", nombre_banco, "] Observaciones diarias:", nrow(rd_xts), "\n")
  list(rd_xts = rd_xts, fechas = as.Date(index(rd_xts)))
}
bbva_diario <- preparar_rend_diario(bbva$r_d, "BBVA")
san_diario  <- preparar_rend_diario(san$r_d,  "SAN")

# --- 11.2 Ajuste y seleccion del mejor GARCH ---
ajustar_garch <- function(diario, nombre_banco) {
  rd <- as.numeric(diario$rd_xts); fechas <- diario$fechas
  cat("\n==================================================\n")
  cat("GARCH para:", nombre_banco, "| Observaciones:", length(rd), "\n")

  # Definimos las cuatro especificaciones (media constante en todas).
  specs <- list(
    sGARCH_norm = ugarchspec(variance.model = list(model = "sGARCH",   garchOrder = c(1, 1)),
                             mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "norm"),
    sGARCH_t    = ugarchspec(variance.model = list(model = "sGARCH",   garchOrder = c(1, 1)),
                             mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std"),
    GJR_GARCH_t = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                             mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std"),
    EGARCH_t    = ugarchspec(variance.model = list(model = "eGARCH",   garchOrder = c(1, 1)),
                             mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
  )

  # Ajustamos cada modelo; si alguno falla, devolvemos NULL en su lugar.
  ajustar_uno <- function(spec_obj) {
    tryCatch(ugarchfit(spec = spec_obj, data = rd, solver = "hybrid", solver.control = list(trace = 0)),
             error = function(e) NULL)
  }
  fits <- lapply(specs, ajustar_uno)
  if (all(sapply(fits, is.null))) stop("Ningun GARCH pudo estimarse para ", nombre_banco)

  # Tabla con AIC, BIC y convergencia de cada modelo.
  ic_tabla <- data.frame(Modelo = names(fits), AIC = NA_real_, BIC = NA_real_,
                         convergencia = NA_integer_, stringsAsFactors = FALSE)
  for (i in seq_along(fits)) if (!is.null(fits[[i]])) {
    ic <- infocriteria(fits[[i]])
    ic_tabla$AIC[i] <- as.numeric(ic[1]); ic_tabla$BIC[i] <- as.numeric(ic[2])
    ic_tabla$convergencia[i] <- fits[[i]]@fit$convergence
  }
  cat("\n--- Comparacion de modelos GARCH ---\n"); print(ic_tabla)

  # Elegimos el mejor por BIC entre los que han convergido (convergencia == 0).
  idx_validos <- which(!is.na(ic_tabla$BIC) & ic_tabla$convergencia == 0)
  if (length(idx_validos) == 0) stop("Ningun GARCH convergio para ", nombre_banco)
  mejor_idx <- idx_validos[which.min(ic_tabla$BIC[idx_validos])]
  mejor_fit <- fits[[mejor_idx]]; mejor_mod <- ic_tabla$Modelo[mejor_idx]
  cat("\nMejor modelo segun BIC:", mejor_mod, "\n")
  print(mejor_fit)

  # Volatilidad condicional estimada y residuos estandarizados z_t.
  vol_cond <- as.numeric(sigma(mejor_fit))
  res_std  <- as.numeric(residuals(mejor_fit, standardize = TRUE))

  # Diagnostico: si el modelo es bueno, z_t y z_t^2 no deben tener autocorrelacion
  # ni efectos ARCH. Como la distribucion es t, no usamos normalidad como criterio;
  # en su lugar comprobamos el ajuste con la transformacion PIT y un test KS.
  lb_res  <- box_ljung_seguro(res_std,   12, 0)
  lb_res2 <- box_ljung_seguro(res_std^2, 12, 0)
  arch_lm <- safe_test(FinTS::ArchTest(res_std, lags = 12))
  lb_res_p <- extraer_p(lb_res); lb_res2_p <- extraer_p(lb_res2); arch_p <- extraer_p(arch_lm)

  shape <- tryCatch(as.numeric(coef(mejor_fit)["shape"]), error = function(e) NA_real_)
  if (grepl("_t$", mejor_mod) && !is.na(shape)) {
    pit <- rugarch::pdist("std", q = res_std, mu = 0, sigma = 1, shape = shape)
    q_teoricos <- rugarch::qdist("std", p = ppoints(length(res_std)), mu = 0, sigma = 1, shape = shape)
  } else {
    pit <- pnorm(res_std); q_teoricos <- qnorm(ppoints(length(res_std)))
  }
  pit <- pmin(pmax(pit, 1e-8), 1 - 1e-8)
  ks_pit_p <- extraer_p(safe_test(stats::ks.test(pit, "punif")))

  cat("\n--- Diagnostico del mejor GARCH ---\n")
  cat("Ljung-Box z_t   lag 12 p-value:", round(lb_res_p,  4), "\n")
  cat("Ljung-Box z_t^2 lag 12 p-value:", round(lb_res2_p, 4), "\n")
  cat("ARCH-LM         lag 12 p-value:", round(arch_p,    4), "\n")
  cat("PIT-KS frente a la distribucion estimada p-value:", round(ks_pit_p, 4), "\n")

  # Grafico: rendimientos y volatilidad condicional estimada.
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  plot(fechas, rd, type = "l", main = paste(nombre_banco, "- Rendimientos diarios"), xlab = "Fecha", ylab = "r_d (%)"); abline(h = 0)
  plot(fechas, vol_cond, type = "l", main = paste(nombre_banco, "- Volatilidad condicional:", mejor_mod), xlab = "Fecha", ylab = "sigma_t (%)")
  par(mfrow = c(1, 1))

  # Grafico de diagnostico de los residuos estandarizados.
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(res_std, type = "l", main = paste(nombre_banco, "- Residuos estandarizados"), xlab = "Observacion", ylab = "z_t"); abline(h = 0); abline(h = c(-3, 3), lty = 2)
  acf(res_std^2, main = paste(nombre_banco, "- ACF de z_t^2"))
  hist(res_std, breaks = 60, probability = TRUE, main = paste(nombre_banco, "- Histograma z_t"), xlab = "z_t")
  if (grepl("_t$", mejor_mod) && !is.na(shape)) curve(rugarch::ddist("std", x, mu = 0, sigma = 1, shape = shape), add = TRUE, lwd = 2)
  else curve(dnorm(x), add = TRUE, lwd = 2)
  qqplot(q_teoricos, sort(res_std), main = paste(nombre_banco, "- QQ-plot vs distribucion estimada"),
         xlab = "Cuantiles teoricos", ylab = "Cuantiles muestrales"); abline(0, 1, lwd = 2)
  par(mfrow = c(1, 1))

  # Prevision de la volatilidad a 20 sesiones.
  prev_garch <- ugarchforecast(mejor_fit, n.ahead = 20)
  plot(1:20, as.numeric(sigma(prev_garch)), type = "l",
       main = paste(nombre_banco, "- Prevision de volatilidad GARCH"), xlab = "Horizonte diario", ylab = "sigma prevista (%)")

  resumen <- data.frame(
    banco = nombre_banco, mejor_modelo_BIC = mejor_mod,
    AIC = round(ic_tabla$AIC[mejor_idx], 5), BIC = round(ic_tabla$BIC[mejor_idx], 5),
    LB_z_p12 = round(lb_res_p, 5), LB_z2_p12 = round(lb_res2_p, 5),
    ARCH_LM_p12 = round(arch_p, 5), PIT_KS_p = round(ks_pit_p, 5), stringsAsFactors = FALSE)

  list(mejor_fit = mejor_fit, mejor_mod = mejor_mod, ic_tabla = ic_tabla,
       vol_cond = vol_cond, res_std = res_std, resumen = resumen,
       vol_df = data.frame(date = fechas, sigma_t = vol_cond))
}

# --- 11.3 Ejecucion ---
cat("\n============================================================\n")
cat("RESULTADOS GARCH\n")
cat("============================================================\n")

garch_bbva <- ajustar_garch(bbva_diario, "BBVA")
garch_san  <- ajustar_garch(san_diario,  "SAN")

tabla_garch <- rbind(garch_bbva$resumen, garch_san$resumen)
rownames(tabla_garch) <- NULL
cat("\nTABLA RESUMEN GARCH:\n"); print(tabla_garch)

# Grafico comparado de la volatilidad condicional de los dos bancos.
n_min     <- min(nrow(garch_bbva$vol_df), nrow(garch_san$vol_df))
bbva_comp <- tail(garch_bbva$vol_df, n_min); san_comp <- tail(garch_san$vol_df, n_min)
plot(bbva_comp$date, bbva_comp$sigma_t, type = "l", main = "Volatilidad condicional GARCH - BBVA vs SAN", xlab = "Fecha", ylab = "sigma_t (%)")
lines(san_comp$date, san_comp$sigma_t, lty = 2)
legend("topright", legend = c("BBVA", "SAN"), lty = c(1, 2), bty = "n")


# --- 11.4 VALIDACION FUERA DE MUESTRA DE LA VOLATILIDAD (rolling con ugarchroll) ---
# Hasta aqui el GARCH se ha elegido y diagnosticado DENTRO de la muestra.
# Este bloque comprueba si el EGARCH-t tiene ademas capacidad PREDICTIVA real:
# reservamos los ultimos dias como test y, avanzando dia a dia, predecimos la
# volatilidad de manana con la informacion disponible hasta hoy. Comparamos esa
# prevision con la volatilidad realmente observada (|r| como proxy) y con dos
# modelos ingenuos de referencia. Si el EGARCH-t tiene menor error que ellos,
# la persistencia y la asimetria detectadas NO son solo ajuste dentro de muestra.
#
# IMPLEMENTACION: usamos ugarchroll, la funcion de rugarch disenada para esto.
# A diferencia de un bucle con ugarchforecast sobre un ajuste fijo (que congelaria
# el origen de prevision entre reestimaciones), ugarchroll hace un rolling 1-dia
# GENUINO: reestima los coeficientes cada 'refit' sesiones y, entre reestimaciones,
# ACTUALIZA el filtrado de la varianza con cada nuevo dato. En cada dia t del tramo
# de test produce sigma_{t} usando solo informacion hasta t-1.
#
# Genera la TABLA 5.11 (tabla_oos) de la memoria.

# Parametros del ejercicio:
#   n_oos : numero de dias reservados como test (los ultimos de la muestra).
#   refit : cada cuantos dias se REESTIMAN los coeficientes (20 ~ 1 mes bursatil).
n_oos_vol <- 500
refit_vol <- 20

# Especificacion EGARCH-t: la MISMA que gana por BIC en el bloque 11.
spec_egarch_t_oos <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

# Funcion que evalua un banco. Devuelve una fila por modelo (EGARCH-t y 2 benchmarks).
evaluar_oos_volatilidad <- function(r_d_xts, nombre_banco,
                                    n_oos = 500, refit = 20) {

  rd <- as.numeric(na.omit(r_d_xts)) * 100   # a %, igual que en el bloque 11
  N  <- length(rd)
  if (N <= n_oos + 300) {
    cat("[", nombre_banco, "] Serie demasiado corta para n_oos =", n_oos, "\n")
    return(NULL)
  }

  cat("\n[", nombre_banco, "] Rolling con ugarchroll:",
      n_oos, "dias de test, refit cada", refit, "(ventana recursiva)...\n")

  # Rolling 1-dia-ahead con reestimacion periodica y filtrado entre reestimaciones.
  # refit.window = "recursive": ventana expansiva (usa todos los datos desde el
  # inicio hasta el punto de reestimacion), coherente con el enfoque del bloque 11.
  roll <- tryCatch(
    ugarchroll(spec_egarch_t_oos, data = rd, n.ahead = 1,
               forecast.length = n_oos, refit.every = refit,
               refit.window = "recursive", solver = "hybrid",
               calculate.VaR = FALSE, keep.coef = FALSE),
    error = function(e) NULL)

  if (is.null(roll)) {
    cat("[", nombre_banco, "] ugarchroll fallo; revisar solver o reducir n_oos.\n")
    return(NULL)
  }

  # Si alguna ventana no convergio, reintentar SOLO esas con otro solver.
  if (!is.null(roll@model$noncidx) && length(roll@model$noncidx) > 0) {
    cat("[", nombre_banco, "] Ventanas sin convergencia:",
        length(roll@model$noncidx), "- reintentando con gosolnp...\n")
    roll <- tryCatch(resume(roll, solver = "gosolnp"), error = function(e) roll)
  }

  # Extraemos las previsiones 1-dia (Sigma) del objeto roll.
  df <- as.data.frame(roll, which = "density")   # columnas: Mu, Sigma, ..., Realized
  sigma_egarch <- as.numeric(df$Sigma)           # prevision 1-dia EGARCH-t (en %)

  # Indices globales del tramo de test (los ultimos n_oos dias) y benchmarks.
  test_idx     <- (N - n_oos + 1):N
  realizado    <- rd[test_idx]                   # rendimiento realizado (en %)
  sigma_rw     <- abs(rd[test_idx - 1])          # benchmark: |r de ayer|
  sigma_uncond <- rep(sd(rd[1:(N - n_oos)]), n_oos)  # benchmark: sd entrenamiento

  # Alineacion robusta por si ugarchroll devolviera menos filas que n_oos.
  m <- min(length(sigma_egarch), n_oos)
  sigma_egarch <- tail(sigma_egarch, m)
  realizado    <- tail(realizado,    m)
  sigma_rw     <- tail(sigma_rw,     m)
  sigma_uncond <- tail(sigma_uncond, m)

  # Proxy de volatilidad realizada: |r| (mismo objetivo para los tres modelos).
  proxy_abs <- abs(realizado)

  mae  <- function(pred, obs) mean(abs(pred - obs), na.rm = TRUE)
  rmse <- function(pred, obs) sqrt(mean((pred - obs)^2, na.rm = TRUE))

  out <- data.frame(
    Banco  = nombre_banco,
    Modelo = c("EGARCH-t (rolling)",
               "Vol. incondicional (constante)",
               "Paseo aleatorio en volatilidad absoluta"),
    MAE  = c(mae(sigma_egarch, proxy_abs),
             mae(sigma_uncond, proxy_abs),
             mae(sigma_rw,     proxy_abs)),
    RMSE = c(rmse(sigma_egarch, proxy_abs),
             rmse(sigma_uncond, proxy_abs),
             rmse(sigma_rw,     proxy_abs)),
    stringsAsFactors = FALSE)
  out$MAE  <- round(out$MAE,  4)
  out$RMSE <- round(out$RMSE, 4)
  out
}

cat("\n============================================================\n")
cat("VALIDACION FUERA DE MUESTRA DE LA VOLATILIDAD (rolling)\n")
cat("============================================================\n")

oos_bbva <- evaluar_oos_volatilidad(bbva$r_d, "BBVA", n_oos = n_oos_vol, refit = refit_vol)
oos_san  <- evaluar_oos_volatilidad(san$r_d,  "SAN",  n_oos = n_oos_vol, refit = refit_vol)

tabla_oos <- rbind(oos_bbva, oos_san)
cat("\nTABLA 5.11 - Evaluacion fuera de muestra de la volatilidad:\n")
print(tabla_oos)

# Lectura automatica: para cada banco, ¿el EGARCH-t tiene el menor MAE y RMSE?
if (!is.null(tabla_oos)) {
  for (b in unique(tabla_oos$Banco)) {
    sub <- tabla_oos[tabla_oos$Banco == b, ]
    gana_mae  <- which.min(sub$MAE)  == 1   # fila 1 = EGARCH-t
    gana_rmse <- which.min(sub$RMSE) == 1
    cat("[", b, "] EGARCH-t bate a los benchmarks -> MAE:", gana_mae,
        "| RMSE:", gana_rmse, "\n")
  }
}


###############################################################
# 12) COINTEGRACION
###############################################################
# Estudiamos relaciones de LARGO PLAZO entre los log-precios bancarios y las
# variables macro en log-niveles:
#   - Mensual:    log(precio banco) con log(IPC)
#   - Trimestral: log(precio banco) con log(PIB)
# Solo tiene sentido si ambas series son I(1). El procedimiento es:
#   1) Comprobar que ambas son I(1) (ADF, PP y KPSS en nivel y en diferencia).
#   2) Engle-Granger: regresion de largo plazo + ADF sobre los residuos.
#   3) Phillips-Ouliaris y Johansen como contrastes especificos de cointegracion.

# --- 12.1 Construccion de los log-niveles ---
# Precios mensuales y trimestrales (ultimo precio de cada periodo).
bbva_px_m <- to.monthly(bbva$px_d, indexAt = "lastof", OHLC = FALSE)
san_px_m  <- to.monthly(san$px_d,  indexAt = "lastof", OHLC = FALSE)
bbva_px_m_df <- data.frame(ym = format(index(bbva_px_m), "%Y-%m"), px_bbva = as.numeric(bbva_px_m))
san_px_m_df  <- data.frame(ym = format(index(san_px_m),  "%Y-%m"), px_san  = as.numeric(san_px_m))

bbva_px_q_raw <- apply.quarterly(bbva$px_d, last)
san_px_q_raw  <- apply.quarterly(san$px_d,  last)
bbva_px_q_df <- data.frame(qkey = qkey_from_date(index(bbva_px_q_raw)), px_bbva = as.numeric(bbva_px_q_raw))
san_px_q_df  <- data.frame(qkey = qkey_from_date(index(san_px_q_raw)),  px_san  = as.numeric(san_px_q_raw))

ipc_m_df <- data.frame(ym   = ym_from_date(ipc_model$Fecha), ipc_val = as.numeric(ipc_model$Valor))
pib_q_df <- data.frame(qkey = qkey_from_date(pib$Fecha),     pib_val = as.numeric(pib$Valor))

# Unimos banco + macro por periodo.
coint_m_bbva <- merge(bbva_px_m_df, ipc_m_df, by = "ym",   all = FALSE)
coint_m_san  <- merge(san_px_m_df,  ipc_m_df, by = "ym",   all = FALSE)
coint_q_bbva <- merge(bbva_px_q_df, pib_q_df, by = "qkey", all = FALSE)
coint_q_san  <- merge(san_px_q_df,  pib_q_df, by = "qkey", all = FALSE)

# Trabajamos en logaritmos (asi los coeficientes se leen como elasticidades).
coint_m_bbva$log_px_bbva <- log(coint_m_bbva$px_bbva); coint_m_bbva$log_ipc <- log(coint_m_bbva$ipc_val)
coint_m_san$log_px_san   <- log(coint_m_san$px_san);   coint_m_san$log_ipc  <- log(coint_m_san$ipc_val)
coint_q_bbva$log_px_bbva <- log(coint_q_bbva$px_bbva); coint_q_bbva$log_pib <- log(coint_q_bbva$pib_val)
coint_q_san$log_px_san   <- log(coint_q_san$px_san);   coint_q_san$log_pib  <- log(coint_q_san$pib_val)

cat("[Cointegracion] Bases de log-niveles construidas.\n")


# --- 12.2 Clasificacion de integracion I(1) ---
# Una serie es I(1) si en nivel NO es estacionaria pero en primera diferencia SI.
# Aplicamos ADF, PP (H0: raiz unitaria) y KPSS (H0: estacionariedad).
#
# IMPORTANTE: para decidir NO usamos una regla que deje a KPSS vetar el resultado.
# Motivo: el p-valor de kpss.test esta truncado en [0.01, 0.1] y pierde potencia
# frente a raices unitarias con deriva, mientras que ADF y PP coinciden con la
# teoria (un log-precio bursatil o un log-indice macro son I(1) de manual). Por eso
# clasificamos por MAYORIA: la serie es I(1) si al menos 2 de los 3 contrastes lo apoyan.
test_integracion <- function(x, nombre) {
  x <- na.omit(as.numeric(x))
  adf_nivel  <- safe_test(tseries::adf.test(x));        adf_diff  <- safe_test(tseries::adf.test(diff(x)))
  pp_nivel   <- safe_test(tseries::pp.test(x));         pp_diff   <- safe_test(tseries::pp.test(diff(x)))
  kpss_nivel <- safe_test(tseries::kpss.test(x, null = "Level")); kpss_diff <- safe_test(tseries::kpss.test(diff(x), null = "Level"))
  adf_nivel_p  <- extraer_p(adf_nivel);  adf_diff_p  <- extraer_p(adf_diff)
  pp_nivel_p   <- extraer_p(pp_nivel);   pp_diff_p   <- extraer_p(pp_diff)
  kpss_nivel_p <- extraer_p(kpss_nivel); kpss_diff_p <- extraer_p(kpss_diff)

  # Cada contraste "vota" I(1) si en nivel apunta a raiz unitaria y en diferencia a estacionariedad.
  voto_adf  <- !is.na(adf_nivel_p)  && !is.na(adf_diff_p)  && adf_nivel_p  > 0.05 && adf_diff_p  < 0.05
  voto_pp   <- !is.na(pp_nivel_p)   && !is.na(pp_diff_p)   && pp_nivel_p   > 0.05 && pp_diff_p   < 0.05
  voto_kpss <- !is.na(kpss_nivel_p) && !is.na(kpss_diff_p) && kpss_nivel_p < 0.05 && kpss_diff_p > 0.05

  # Clasificacion por mayoria (>= 2 de 3). KPSS cuenta como un voto, no como un veto.
  votos_i1 <- sum(voto_adf, voto_pp, voto_kpss, na.rm = TRUE)
  es_i1    <- votos_i1 >= 2

  cat("\n--------------------------------------------------\n")
  cat("Integracion:", nombre, "\n")
  cat("ADF  nivel p:", round(adf_nivel_p,  4), " | ADF  diff p:", round(adf_diff_p,  4), "\n")
  cat("PP   nivel p:", round(pp_nivel_p,   4), " | PP   diff p:", round(pp_diff_p,   4), "\n")
  cat("KPSS nivel p:", round(kpss_nivel_p, 4), " | KPSS diff p:", round(kpss_diff_p, 4), "\n")
  cat("Votos a favor de I(1):", votos_i1, "de 3  ->  Clasificacion:", ifelse(es_i1, "I(1)", "No I(1)"), "\n")

  data.frame(serie = nombre,
             ADF_nivel_p = round(adf_nivel_p, 5), ADF_diff_p = round(adf_diff_p, 5),
             PP_nivel_p  = round(pp_nivel_p,  5), PP_diff_p  = round(pp_diff_p,  5),
             KPSS_nivel_p = round(kpss_nivel_p, 5), KPSS_diff_p = round(kpss_diff_p, 5),
             votos_I1 = votos_i1, es_I1 = es_i1, stringsAsFactors = FALSE)
}


# --- 12.3 Engle-Granger ---
# Paso 1: regresion de largo plazo  y = b0 + b1*x + u.
# Paso 2: ADF sobre los residuos u (sin constante, porque tienen media cero).
# Este ADF se usa solo como diagnostico auxiliar: la decision formal de
# cointegracion se apoya en Phillips-Ouliaris y Johansen.
test_engle_granger <- function(y, x, nombre_y, nombre_x) {
  y <- na.omit(as.numeric(y)); x <- na.omit(as.numeric(x))
  n <- min(length(y), length(x)); y <- tail(y, n); x <- tail(x, n)
  cat("\n--------------------------------------------------\n")
  cat("Engle-Granger:", nombre_y, "~", nombre_x, "\n")
  reg_lp <- lm(y ~ x); print(summary(reg_lp))
  res_lp <- residuals(reg_lp)
  adf_res <- urca::ur.df(res_lp, type = "none", selectlags = "AIC")
  cat("\nADF sobre los residuos de largo plazo:\n"); print(summary(adf_res))
  stat <- as.numeric(adf_res@teststat[1])
  crit5 <- as.numeric(adf_res@cval[1, "5pct"])
  cat("\nEstadistico ADF residuos:", round(stat, 4), "| critico 5% estandar (no decisorio):", round(crit5, 4), "\n")
  list(reg = reg_lp, residuos = res_lp, stat = stat, crit_5 = crit5, y = y, x = x)
}

# --- 12.4 Phillips-Ouliaris ---
# Contraste especifico de cointegracion (rechaza H0 de NO cointegracion si el
# estadistico supera el valor critico al 5%).
test_phillips_ouliaris <- function(y, x, nombre_par) {
  mat <- na.omit(cbind(y = as.numeric(y), x = as.numeric(x)))
  cat("\n--------------------------------------------------\n")
  cat("Phillips-Ouliaris:", nombre_par, "\n")
  po <- urca::ca.po(mat, demean = "constant", type = "Pz"); print(summary(po))
  po
}

# --- 12.5 Johansen ---
# Contraste multivariante (Johansen, 1988) con los tests de la traza y del maximo
# valor propio. Para 2 variables, si se rechaza r = 0 hay al menos un vector de cointegracion.
test_johansen <- function(mat, nombre_par, max_lag = 6) {
  mat <- na.omit(as.matrix(mat))
  cat("\n--------------------------------------------------\n")
  cat("Johansen:", nombre_par, "\n")
  if (nrow(mat) < 30) cat("Muestra reducida; interpretar con cautela.\n")
  # Retardos: como maximo nrow/5, y al menos 2 (lo exige ca.jo).
  lag_max_real <- max(min(max_lag, floor(nrow(mat) / 5)), 2)
  sel <- vars::VARselect(mat, lag.max = lag_max_real, type = "const")
  k   <- as.integer(sel$selection["AIC(n)"]); if (is.na(k) || k < 2) k <- 2
  cat("Retardo K usado:", k, "\n")
  jo_trace <- urca::ca.jo(mat, type = "trace", ecdet = "const", K = k)
  jo_eigen <- urca::ca.jo(mat, type = "eigen", ecdet = "const", K = k)
  cat("\nTest de la traza:\n");             print(summary(jo_trace))
  cat("\nTest del maximo valor propio:\n"); print(summary(jo_eigen))
  list(trace = jo_trace, eigen = jo_eigen, K = k)
}

# Extrae la decision de Johansen para la hipotesis r = 0 (existencia de cointegracion).
extraer_decision_johansen <- function(jo_obj) {
  if (is.null(jo_obj)) {
    return(data.frame(Johansen_trace_coint = NA, Johansen_eigen_coint = NA))
  }
  decision_r0 <- function(obj) {
    rn <- rownames(obj@cval); idx <- grep("r = 0|r=0", rn)
    if (length(idx) == 0) idx <- length(obj@teststat)
    as.numeric(obj@teststat[idx]) > as.numeric(obj@cval[idx, "5pct"])
  }
  data.frame(Johansen_trace_coint = decision_r0(jo_obj$trace),
             Johansen_eigen_coint = decision_r0(jo_obj$eigen))
}

# --- 12.6 Graficos de cointegracion ---
graficar_cointegracion <- function(y, x, res_lp, fecha_vec, nombre_y, nombre_x, frecuencia) {
  n <- min(length(y), length(x), length(res_lp), length(fecha_vec))
  y <- tail(as.numeric(y), n); x <- tail(as.numeric(x), n)
  res_lp <- tail(as.numeric(res_lp), n); fecha_vec <- tail(fecha_vec, n)
  fechas <- if (frecuencia == "mensual") as.Date(paste0(fecha_vec, "-01")) else qkey_to_date(fecha_vec)
  # Normalizamos para comparar ambas series en la misma escala.
  y_norm <- (y - mean(y)) / sd(y); x_norm <- (x - mean(x)) / sd(x)
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
  plot(fechas, y_norm, type = "l", main = paste("Series normalizadas:", nombre_y, "vs", nombre_x), xlab = "Fecha", ylab = "valor normalizado")
  lines(fechas, x_norm, lty = 2); legend("topleft", legend = c(nombre_y, nombre_x), lty = c(1, 2), bty = "n")
  plot(fechas, res_lp, type = "l", main = "Residuo de largo plazo", xlab = "Fecha", ylab = "residuo"); abline(h = 0, lty = 2)
  acf(res_lp, main = "ACF del residuo de largo plazo")
  par(mfrow = c(1, 1))
}

# --- 12.7 Funcion que analiza un par completo ---
analizar_cointegracion_par <- function(df, y_col, x_col, fecha_col,
                                       nombre_y, nombre_x, nombre_par, frecuencia, max_lag_johansen) {
  y <- df[[y_col]]; x <- df[[x_col]]; fecha_vec <- df[[fecha_col]]
  cat("\n==================================================\n")
  cat("ANALISIS DE COINTEGRACION:", nombre_par, "\n")

  # Paso 1: comprobar que ambas series son I(1).
  int_y <- test_integracion(y, nombre_y)
  int_x <- test_integracion(x, nombre_x)
  ambas_i1 <- isTRUE(int_y$es_I1[1]) && isTRUE(int_x$es_I1[1])

  eg <- NULL; po <- NULL; jo <- NULL
  if (ambas_i1) {
    # Paso 2 y 3: contrastes de cointegracion.
    eg <- test_engle_granger(y, x, nombre_y, nombre_x)
    po <- test_phillips_ouliaris(y, x, nombre_par)
    jo <- test_johansen(cbind(y, x), nombre_par, max_lag = max_lag_johansen)
    graficar_cointegracion(y, x, eg$residuos, fecha_vec, nombre_y, nombre_x, frecuencia)
  } else {
    cat("\nNo se aplica cointegracion: al menos una serie no es I(1).\n")
  }

  # Decision de Phillips-Ouliaris.
  po_stat <- if (!is.null(po)) as.numeric(po@teststat) else NA_real_
  po_crit <- if (!is.null(po)) as.numeric(po@cval[1, "5pct"]) else NA_real_
  po_coint <- if (!is.null(po)) po_stat > po_crit else NA
  jo_dec <- extraer_decision_johansen(jo)

  # Contamos en cuantos contrastes hay cointegracion (PO, Johansen traza, Johansen valor propio).
  evidencias <- c(isTRUE(po_coint), isTRUE(jo_dec$Johansen_trace_coint[1]), isTRUE(jo_dec$Johansen_eigen_coint[1]))
  num_positivos  <- if (ambas_i1) sum(evidencias, na.rm = TRUE) else NA_integer_
  coint_algun    <- if (ambas_i1) num_positivos >= 1 else NA   # evidencia en al menos un contraste
  coint_coherente<- if (ambas_i1) num_positivos >= 2 else NA   # evidencia coherente (>= 2 contrastes)

  resumen <- data.frame(
    par = nombre_par, frecuencia = frecuencia, y = nombre_y, x = nombre_x,
    y_I1 = int_y$es_I1[1], x_I1 = int_x$es_I1[1], ambas_I1 = ambas_i1,
    EG_ADF_resid = ifelse(is.null(eg), NA, round(eg$stat, 5)),
    PO_stat = ifelse(ambas_i1, round(po_stat, 5), NA),
    PO_crit5 = ifelse(ambas_i1, round(po_crit, 5), NA),
    PO_cointegracion = ifelse(ambas_i1, po_coint, NA),
    jo_dec,
    N_contrastes_positivos = num_positivos,
    Cointegracion_algun_contraste = coint_algun,
    Cointegracion_coherente = coint_coherente,
    stringsAsFactors = FALSE)

  list(integracion_y = int_y, integracion_x = int_x, resumen = resumen)
}

# --- 12.8 Ejecucion ---
cat("\n============================================================\n")
cat("RESULTADOS COINTEGRACION\n")
cat("============================================================\n")

coint_bbva_ipc <- analizar_cointegracion_par(coint_m_bbva, "log_px_bbva", "log_ipc", "ym",
  "log(BBVA)", "log(IPC)", "BBVA + IPC", "mensual", 12)
coint_san_ipc  <- analizar_cointegracion_par(coint_m_san,  "log_px_san",  "log_ipc", "ym",
  "log(SAN)", "log(IPC)", "SAN + IPC", "mensual", 12)
coint_bbva_pib <- analizar_cointegracion_par(coint_q_bbva, "log_px_bbva", "log_pib", "qkey",
  "log(BBVA)", "log(PIB)", "BBVA + PIB", "trimestral", 6)
coint_san_pib  <- analizar_cointegracion_par(coint_q_san,  "log_px_san",  "log_pib", "qkey",
  "log(SAN)", "log(PIB)", "SAN + PIB", "trimestral", 6)

tabla_coint <- rbind(coint_bbva_ipc$resumen, coint_san_ipc$resumen,
                     coint_bbva_pib$resumen, coint_san_pib$resumen)
tabla_integracion_coint <- rbind(
  coint_bbva_ipc$integracion_y, coint_bbva_ipc$integracion_x,
  coint_san_ipc$integracion_y,  coint_san_ipc$integracion_x,
  coint_bbva_pib$integracion_y, coint_bbva_pib$integracion_x,
  coint_san_pib$integracion_y,  coint_san_pib$integracion_x
)

cat("\nTABLA RESUMEN COINTEGRACION:\n"); print(tabla_coint)
cat("\nTABLA DE INTEGRACION (LOG-NIVELES):\n"); print(tabla_integracion_coint)


###############################################################
# 13) RESUMEN FINAL
###############################################################

cat("\n============================================================\n")
cat("RESUMEN FINAL DE MODELOS\n")
cat("============================================================\n")
cat("\n1) ARIMA / SARIMA\n"); print(tabla_arima)
cat("\n2) GARCH\n");          print(tabla_garch)
cat("\n3) COINTEGRACION\n");  print(tabla_coint)
cat("\n4) GRANGER\n");        print(tabla_granger)
if (exists("tabla_oos") && !is.null(tabla_oos)) {
  cat("\n5) VALIDACION FUERA DE MUESTRA (volatilidad)\n"); print(tabla_oos)
}

# Guardamos las tablas principales en CSV por si las necesitamos para la memoria.
write.csv(tabla_arima,             "tabla_arima.csv",                  row.names = FALSE)
write.csv(tabla_garch,             "tabla_garch.csv",                  row.names = FALSE)
write.csv(tabla_coint,             "tabla_cointegracion.csv",          row.names = FALSE)
write.csv(tabla_integracion_coint, "tabla_integracion_cointegracion.csv", row.names = FALSE)
write.csv(tabla_granger,           "tabla_granger.csv",                row.names = FALSE)
write.csv(tabla_granger_sensibilidad_ipc, "tabla_granger_sensibilidad_ipc.csv", row.names = FALSE)
write.csv(tabla_granger_sensibilidad_pib, "tabla_granger_sensibilidad_pib.csv", row.names = FALSE)
write.csv(tabla_granger_pib_pre2020,      "tabla_granger_pib_pre2020.csv",      row.names = FALSE)
write.csv(tabla_crisis,            "tabla_crisis.csv",                 row.names = FALSE)
if (exists("tabla_oos") && !is.null(tabla_oos)) {
  write.csv(tabla_oos,             "tabla_oos_volatilidad.csv",        row.names = FALSE)
}
cat("\nTablas guardadas en CSV.\n")

# Cerramos el PDF con todas las figuras.
if (dev.cur() > 1) dev.off()

# Informacion de la sesion (versiones de R y de los paquetes) para reproducibilidad.
# Conviene ejecutar el script en una sesion limpia y conservar esta salida.
cat("\n============================================================\n")
cat("INFORMACION DE LA SESION (sessionInfo)\n")
cat("============================================================\n")
print(sessionInfo())
