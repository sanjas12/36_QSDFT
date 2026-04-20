#!/bin/bash

# uv + интернет → uv sync (если есть pyproject.toml) или uv pip install -r requirements.txt
# uv + офлайн → uv sync --no-index --find-links=... или uv pip install --no-index
# нет uv + интернет → fallback на pip install -r requirements.txt
# нет uv + офлайн → fallback на pip install --no-index -f ...


# Переходим в корень проекта (на уровень выше от scripts/)
cd "$(dirname "$0")/.." || exit 1

LOCAL_PACKAGES_DIR="d:\\temp\\python_Library"


check_internet() {
    local result=0

    if command -v uv &>/dev/null; then
        echo "🌐 Проверка интернета через uv..."
        #  --refresh принудительно обращается к PyPI, игнорируя кэш
        uv pip install pip --dry-run --refresh --quiet 2>/dev/null
        result=$?
    elif command -v curl &>/dev/null; then
        echo "🌐 Проверка интернета через curl → pypi.org..."
        curl -s --connect-timeout 3 --max-time 5 https://pypi.org >/dev/null 2>&1
        result=$?
    elif command -v wget &>/dev/null; then
        echo "🌐 Проверка интернета через wget → pypi.org..."
        wget -q --timeout=3 --tries=1 https://pypi.org -O /dev/null 2>&1
        result=$?
    else
        echo "⚠️  Нет доступных инструментов для проверки (uv/curl/wget)"
        return 1
    fi

    if [ $result -eq 0 ]; then
        echo "✅ Интернет доступен"
    else
        echo "❌ Интернет недоступен (код: $result)"
    fi

    return $result
}


activate_venv() {
    if [[ -f ".venv/Scripts/activate" ]]; then
        source .venv/Scripts/activate   # Windows Git Bash
    elif [[ -f ".venv/bin/activate" ]]; then
        source .venv/bin/activate       # Linux / macOS
    else
        echo "❌ Не удалось найти activate в .venv"
        exit 1
    fi
}

# ── Проверка наличия uv ───────────────────────────────────────────────────────

if ! command -v uv &>/dev/null; then
    echo "⚠️  uv не найден — используем fallback на pip"

    if ! check_internet; then
        echo "❌ Интернет недоступен и uv не установлен — локальная установка через pip"
        rm -rf .venv
        python -m venv .venv
        activate_venv

        if [ ! -f "requirements.txt" ]; then
            echo "❌ Файл requirements.txt не найден"
            exit 1
        fi

        pip install -r requirements.txt --no-index -f "$LOCAL_PACKAGES_DIR" --no-deps
        echo "✅ .venv создан через pip (локально, без uv)"
    else
        echo "Интернет доступен — устанавливаем через pip"
        rm -rf .venv
        python -m venv .venv
        activate_venv

        if [ ! -f "requirements.txt" ]; then
            echo "❌ Файл requirements.txt не найден"
            exit 1
        fi

        python -m pip install --upgrade pip
        pip install -r requirements.txt
        echo "✅ .venv создан через pip (uv не найден)"
    fi

    exit 0
fi

# ── uv найден ─────────────────────────────────────────────────────────────────

echo "✅ uv найден: $(uv --version)"

# ── Офлайн-режим ─────────────────────────────────────────────────────────────
if ! check_internet; then
    echo "❌ PyPi.org недоступен — офлайн-установка через uv"

    if [ ! -d "$LOCAL_PACKAGES_DIR" ] && [[ "$LOCAL_PACKAGES_DIR" != /* ]]; then
        # Windows-путь передан — uv понимает его только в нативном виде
        UV_CACHE_FLAG="--find-links=$(cygpath -u "$LOCAL_PACKAGES_DIR" 2>/dev/null || echo "$LOCAL_PACKAGES_DIR")"
    else
        UV_CACHE_FLAG="--find-links=$LOCAL_PACKAGES_DIR"
    fi

    rm -rf .venv

    if [ -f "pyproject.toml" ]; then
        uv sync --no-index $UV_CACHE_FLAG
    elif [ -f "requirements.txt" ]; then
        uv venv
        uv pip install -r requirements.txt --no-index $UV_CACHE_FLAG
    else
        echo "❌ Не найден ни pyproject.toml, ни requirements.txt"
        exit 1
    fi

    activate_venv
    echo "✅ .venv создан через uv (офлайн)"
    exit 0
fi

# ── Онлайн-режим: uv есть, интернет есть ─────────────────────────────────────
rm -rf .venv

if [ -f "pyproject.toml" ]; then
    echo "Найден pyproject.toml — используем uv sync"
    uv sync
elif [ -f "requirements.txt" ]; then
    echo "Найден requirements.txt — используем uv pip install"
    uv venv
    uv pip install -r requirements.txt
else
    echo "❌ Не найден ни pyproject.toml, ни requirements.txt"
    exit 1
fi

activate_venv
echo "✅ .venv создан и активирован через uv"