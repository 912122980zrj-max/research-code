@echo off
chcp 65001 >nul
cd /d "%~dp0"

set "RSCRIPT=C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
if not exist "%RSCRIPT%" (
  where Rscript >nul 2>nul
  if errorlevel 1 (
    echo [ERROR] Rscript 未找到，请先安装 R 或修改本脚本中的 RSCRIPT 路径。
    pause
    exit /b 1
  )
  set "RSCRIPT=Rscript"
)

echo ============================================
echo  开始一键安装 TBI 缺失 R 包（可安全重复运行）
echo ============================================
"%RSCRIPT%" --vanilla "%~dp0install_missing_packages.R"

echo.
echo 安装完成，请查看 INSTALL_REPORT.txt 和 INSTALL_LOG_*.txt
pause
