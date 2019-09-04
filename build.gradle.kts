import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    java
    kotlin("jvm") version "1.3.41"
    id("org.jetbrains.kotlin.plugin.noarg") version "1.3.41"
    id("org.jetbrains.kotlin.plugin.allopen") version "1.3.41"
}

group = "io.github.shinglem"

repositories {
    mavenCentral()
}

dependencies {
    implementation(kotlin("stdlib-jdk8"))
    testCompile("junit", "junit", "4.12")
}

configure<JavaPluginConvention> {
    sourceCompatibility = JavaVersion.VERSION_1_8
}
tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}

subprojects {

    repositories {
        mavenLocal()
        maven("http://maven.aliyun.com/nexus/content/groups/public/")
        maven("https://dl.bintray.com/kotlin/exposed")
        maven("https://jitpack.io")
        jcenter()
        mavenCentral()
    }
}