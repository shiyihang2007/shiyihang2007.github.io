---
title: vscode-tasks
date: 2023-11-17 10:19:21
---

善用 Tasks 可以救命喵！

### 考场用简化版本

#### 测试程序运行时间

```json
{
	"version": "2.0.0",
	"tasks": [
		{
			"type": "process",
			"label": "C/C++ (g++, O2) 生成活动文件",
			"command": "g++",
			"args": [
				"-fdiagnostics-color=always",
				"${file}",
				"-std=c++17",
				"-Wall",
				"-Wextra",
				"-O2",
				"-DTEST",
				"-o",
				"${fileDirname}/${fileBasenameNoExtension}"
			],
			"problemMatcher": "$gcc"
		},
		{
			"label": "Test Usage of Program",
			"type": "process",
			"command": "/usr/bin/time",
			"args": [
				"-f========\\nTime Usage\\n%es real\\n%Us user\\n%Ss sys\\nMemory Usage\\n%MKB set\\n========",
				"${fileDirname}/${fileBasenameNoExtension}"
			],
			"options": {
				"cwd": "${fileDirname}"
			},
			"dependsOn": ["C/C++ (g++, O2) 生成活动文件"]
		}
	]
}

```

#### 调试

注意： **需要 C/C++ 插件支持**

tasks.json
```json
{
	"version": "2.0.0",
	"tasks": [
		{
			"type": "process",
			"label": "C/C++ (g++, O2) 生成活动文件",
			"command": "g++",
			"args": [
				"-fdiagnostics-color=always",
				"${file}",
				"-std=c++17",
				"-Wall",
				"-Wextra",
				"-O2",
				"-DTEST",
				"-o",
				"${fileDirname}/${fileBasenameNoExtension}"
			],
			"problemMatcher": "$gcc"
		}
	]
}

```

launch.json
```json
{
	"configurations": [
		{
			"name": "C/C++ (gdb) 启动",
			"type": "cppdbg",
			"request": "launch",
			"args": [],
			"stopAtEntry": false,
			"cwd": "${fileDirname}",
			"externalConsole": false,
			"miDebuggerPath": "/usr/bin/gdb",
			"program": "${fileDirname}/${fileBasenameNoExtension}",
			"preLaunchTask": "C/C++ (g++) 生成活动文件",
			"MIMode": "gdb"
		}
	]
}

```

### 个人自用版本

```json
// tasks.json, 放在 用户数据文件夹 或 项目根目录/.vscode/ 下
{
	"version": "2.0.0",  // 解析器版本号，不要改
	"tasks": [           // 任务列表
		{                // 这里的每一项都会作为一个任务显示
			"type": "process",                   // 作为进程(process)或命令行(shell)处理
			"label": "C/C++ (g++) 生成活动文件",  // 任务名称
			"command": "g++",                    // 要运行的程序
			"linux": {                           // 对linux系统特化选项
				"args": [                        // 程序运行参数
					"-fdiagnostics-color=always",
					"-g",
					"${file}",
					"-std=c++17",
					"-Wall",
					"-Wextra",
					"-DTEST",
					"-o",
					"${fileDirname}/${fileBasenameNoExtension}"
				],
			},
			"windows": {                        // 对windows系统特化选项
				"args": [                       // 程序运行参数
					"-fdiagnostics-color=always",
					"-g",
					"${file}",
					"-std=c++17",
					"-lm",
					"-Wall",
					"-Wextra",
					"-DTEST",
					"-o",
					"${fileDirname}\\${fileBasenameNoExtension}.exe"
				],
			},
			"problemMatcher": {                 // 匹配错误
				"owner": "cpp",                 // 仅对于 C++ 文件
				"fileLocation": [               // 匹配文件的位置
					"relative",                 // 相对位置
					"${workspaceRoot}"          // 相对于工作区根目录
				],
				"pattern": {                    // 匹配方式
					"regexp": "^(.*):(\\d+):(\\d+):\\s+(warning|error):\\s+(.*)$",
						// 自己对着 g++ 报错理解一下，这里的括号和下面的序号对应
					"file": 1,     // 第一个括号匹配文件
					"line": 2,     // 行号
					"column": 3,   // 列号
					"severity": 4, // 严重程度
					"message": 5   // 详细信息
				}
			},
			"group": {                          // 作为什么任务类型
				"kind": "build",                // 可以是编译(build)和测试(test)
				"isDefault": true               // 是否作为默认任务运行
			},
			"icon": {                           // 任务的图标，会显示为内置终端的图标
				"id": "window",                 // 图标样式
				"color": "terminal.ansiWhite"   // 图标颜色
			}
		},
		{
			"type": "process",                  // 同上，但是开 O2
			"label": "C/C++ (g++, O2) 生成活动文件",
			"command": "g++",
			"linux": {
				"args": [
					"-fdiagnostics-color=always",
					"${file}",
					"-std=c++17",
					"-Wall",
					"-Wextra",
					"-O2",
					"-DTEST",
					"-o",
					"${fileDirname}/${fileBasenameNoExtension}"
				],
			},
			"windows": {
				"args": [
					"-fdiagnostics-color=always",
					"${file}",
					"-std=c++17",
					"-lm",
					"-O2",
					"-Wall",
					"-Wextra",
					"-DTEST",
					"-o",
					"${fileDirname}\\${fileBasenameNoExtension}.exe"
				],
			},
			"problemMatcher": {
				"owner": "cpp",
				"fileLocation": [
					"relative",
					"${workspaceRoot}"
				],
				"pattern": {
					"regexp": "^(.*):(\\d+):(\\d+):\\s+(warning|error):\\s+(.*)$",
					"file": 1,
					"line": 2,
					"column": 3,
					"severity": 4,
					"message": 5
				}
			},
			"group": {
				"kind": "build",
				"isDefault": true
			},
			"icon": {
				"id": "window",
				"color": "terminal.ansiWhite"
			}
		},
		{
			"type": "process",                  // 这个任务用来测试文件运行速度
			"label": "Test Usage of Program",
			"command": "/usr/bin/time",
			"args": [
				"-f========\\nTime Usage\\n%es real\\n%Us user\\n%Ss sys\\nMemory Usage\\n%MKB set\\n========",
				"${fileDirname}/${fileBasenameNoExtension}"
			],
			"options": {
				"cwd": "${fileDirname}",
			},
			"dependsOn": [                      // 在运行此任务前会先运行这里的任务，注意这里的任务不指定是乱序执行的
				"C/C++ (g++, O2) 生成活动文件"
			],
			"group": {
				"kind": "test",
				"isDefault": true
			},
			"icon": {
				"id": "dashboard",
				"color": "terminal.ansiWhite"
			}
		}
	]
}

```


